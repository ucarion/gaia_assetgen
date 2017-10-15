#[macro_use]
extern crate error_chain;
#[macro_use]
extern crate serde_derive;

extern crate colored;
extern crate geo;
extern crate geojson;
extern crate serde;
extern crate serde_json;

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{Read, Write};
use std::path::PathBuf;
use std::process::{Command, Stdio};
use std::str::FromStr;

use colored::Colorize;
use geojson::GeoJson;
use geojson::conversion::TryInto;
use geo::boundingbox::BoundingBox;
use geo::simplify::Simplify;

mod errors {
    error_chain!{}
}

use errors::*;

pub const ELEVATION_OFFSET: u16 = 500;

// This is just for logging purposes
const NUM_STEPS: u32 = 6;

const MAX_LEVEL: u32 = 6;
const ELEVATION_TILE_SIZE: u32 = 64;
const PIXELS_PER_TILE: u32 = 256;

const TILES_ACROSS_WIDTH: u32 = 128;
const NASA_TILES_ACROSS_WIDTH: u32 = 4;

// Constants relating to NOAA GLOBE data
const NOAA_TILE_WIDTH: u32 = 10800;
const NOAA_TILE_HEIGHT_TOP: u32 = 4800;
const NOAA_TILE_HEIGHT_MIDDLE: u32 = 6000;

// The metadata for the tile metadata on the first pass; only maxima and minima, computed by
// Imagemagick, are present.
#[derive(Serialize, Deserialize, Debug)]
struct FirstPassTileMetadata {
    min_elevation: u16,
    max_elevation: u16,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct TileMetadata {
    pub min_elevation: u16,
    pub max_elevation: u16,
    pub polygons: Vec<u64>,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct PolygonPointData {
    pub polygons: Vec<MultiLevelPolygon>,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct MultiLevelPolygon {
    pub properties: serde_json::Map<String, serde_json::Value>,
    pub levels: Vec<Vec<(f32, f32)>>,
}

#[derive(Debug)]
pub struct PrepareAssetsTask {
    noaa_globe_dir: PathBuf,
    nasa_blue_marble_dir: PathBuf,
    features_file: PathBuf,
    simplification_epsilons: [f32; MAX_LEVEL as usize + 1],
    output_dir: PathBuf,
}

impl PrepareAssetsTask {
    pub fn new() -> PrepareAssetsTask {
        PrepareAssetsTask {
            noaa_globe_dir: "".into(),
            nasa_blue_marble_dir: "".into(),
            features_file: "".into(),
            simplification_epsilons: [0.0; MAX_LEVEL as usize + 1],
            output_dir: "".into(),
        }
    }

    pub fn with_noaa_globe_dir(self, noaa_globe_dir: PathBuf) -> PrepareAssetsTask {
        PrepareAssetsTask {
            noaa_globe_dir: noaa_globe_dir,
            ..self
        }
    }

    pub fn with_nasa_blue_marble_dir(self, nasa_blue_marble_dir: PathBuf) -> PrepareAssetsTask {
        PrepareAssetsTask {
            nasa_blue_marble_dir: nasa_blue_marble_dir,
            ..self
        }
    }

    pub fn with_features_file(self, features_file: PathBuf) -> PrepareAssetsTask {
        PrepareAssetsTask {
            features_file: features_file,
            ..self
        }
    }

    pub fn with_output_dir(self, output_dir: PathBuf) -> PrepareAssetsTask {
        PrepareAssetsTask {
            output_dir: output_dir,
            ..self
        }
    }

    pub fn with_simplification_epsilons(
        self,
        simplification_epsilons: [f32; MAX_LEVEL as usize + 1],
    ) -> PrepareAssetsTask {
        PrepareAssetsTask {
            simplification_epsilons: simplification_epsilons,
            ..self
        }
    }

    pub fn run(&self) -> Result<()> {
        if !self.crops_dir().is_dir() {
            std::fs::create_dir(&self.crops_dir()).chain_err(
                || "Could not create crops directory",
            )?;
        }

        if !self.tiles_dir().is_dir() {
            std::fs::create_dir(&self.tiles_dir()).chain_err(
                || "Could not create tiles directory",
            )?;
        }

        self.create_nasa_level0_crops()?;
        self.create_noaa_level0_crops()?;

        for level in 1..MAX_LEVEL + 1 {
            self.create_crop_level(level)?;
        }

        for level in 0..MAX_LEVEL + 1 {
            self.create_tile_level(level)?;
        }

        self.create_polygon_data()?;

        Ok(())
    }

    fn create_nasa_level0_crops(&self) -> Result<()> {
        if self.crops_dir().join("0_0_0.jpg").is_file() {
            return Ok(());
        }

        for y in 0..2 {
            for x in 0..4 {
                self.crop_blue_marble_tile(x, y)?;
            }
        }

        Ok(())
    }

    fn crop_blue_marble_tile(&self, x: u32, y: u32) -> Result<()> {
        self.pretty_log(
            1,
            1 + x + y * NASA_TILES_ACROSS_WIDTH,
            NASA_TILES_ACROSS_WIDTH * NASA_TILES_ACROSS_WIDTH / 2,
            &format!("Resizing and cropping Blue Marble tile {}.{}", x, y),
        );

        let name_offset_x = x * TILES_ACROSS_WIDTH / NASA_TILES_ACROSS_WIDTH;
        let name_offset_y = y * TILES_ACROSS_WIDTH / NASA_TILES_ACROSS_WIDTH;
        let name_template_x = format!("%[fx:page.x / {} + {}]", PIXELS_PER_TILE, name_offset_x);
        let name_template_y = format!("%[fx:page.y / {} + {}]", PIXELS_PER_TILE, name_offset_y);
        let name_template = format!("0_{}_{}.jpg", name_template_x, name_template_y);

        let tile_name = format!(
            "{}{}",
            ['A', 'B', 'C', 'D'][x as usize],
            ['1', '2'][y as usize]
        );
        let tile_file_name = format!("world.topo.bathy.200412.3x21600x21600.{}.jpg", tile_name);
        let tile_path = self.nasa_blue_marble_dir.join(tile_file_name);

        let resize_size = PIXELS_PER_TILE * TILES_ACROSS_WIDTH / NASA_TILES_ACROSS_WIDTH;
        let crop_size = format!("{}x{}", PIXELS_PER_TILE, PIXELS_PER_TILE);

        let output = Command::new("convert")
            .arg("-monitor")
            .arg(tile_path)
            .args(&["-resize", &resize_size.to_string()])
            .args(&["-crop", &crop_size])
            .args(&["-set", "filename:tile", &name_template])
            .args(&["+repage", "+adjoin"])
            .arg(self.crops_dir().join("%[filename:tile]"))
            .stderr(Stdio::inherit())
            .output()
            .chain_err(|| "Error while running `convert`")?;

        if !output.status.success() {
            return Err("`convert` returned with non-zero exit status".into());
        }

        Ok(())
    }

    fn create_noaa_level0_crops(&self) -> Result<()> {
        if self.crops_dir().join("0_0_0.elevation").is_file() {
            return Ok(());
        }

        self.pretty_log(2, 1, 4, "Cropping GLOBE tile quadrant 1");
        self.crop_noaa_tile_quadrant(
            (0, 0),
            (NOAA_TILE_HEIGHT_TOP, NOAA_TILE_HEIGHT_MIDDLE),
            ["a10g", "b10g", "e10g", "f10g"],
        )?;
        self.pretty_log(2, 2, 4, "Cropping GLOBE tile quadrant 2");
        self.crop_noaa_tile_quadrant(
            (64, 0),
            (NOAA_TILE_HEIGHT_TOP, NOAA_TILE_HEIGHT_MIDDLE),
            ["c10g", "d10g", "g10g", "h10g"],
        )?;
        self.pretty_log(2, 3, 4, "Cropping GLOBE tile quadrant 3");
        self.crop_noaa_tile_quadrant(
            (0, 32),
            (NOAA_TILE_HEIGHT_MIDDLE, NOAA_TILE_HEIGHT_TOP),
            ["i10g", "j10g", "m10g", "n10g"],
        )?;
        self.pretty_log(2, 4, 4, "Cropping GLOBE tile quadrant 4");
        self.crop_noaa_tile_quadrant(
            (64, 32),
            (NOAA_TILE_HEIGHT_MIDDLE, NOAA_TILE_HEIGHT_TOP),
            ["k10g", "l10g", "o10g", "p10g"],
        )?;
        Ok(())
    }

    fn crop_noaa_tile_quadrant(
        &self,
        name_offset: (u32, u32),
        tile_heights: (u32, u32),
        tiles: [&str; 4],
    ) -> Result<()> {
        let name_template_x = format!("%[fx:page.x / {} + {}]", ELEVATION_TILE_SIZE, name_offset.0);
        let name_template_y = format!("%[fx:page.y / {} + {}]", ELEVATION_TILE_SIZE, name_offset.1);
        let name_template = format!("0_{}_{}.elevation", name_template_x, name_template_y);

        let top_left = self.noaa_globe_dir.join(tiles[0]);
        let top_right = self.noaa_globe_dir.join(tiles[1]);
        let top_size = format!("{}x{}", NOAA_TILE_WIDTH, tile_heights.0);

        let bottom_left = self.noaa_globe_dir.join(tiles[2]);
        let bottom_right = self.noaa_globe_dir.join(tiles[3]);
        let bottom_size = format!("{}x{}", NOAA_TILE_WIDTH, tile_heights.1);

        let offset_resize_and_append = &[
            "-evaluate",
            "addmodulus",
            &ELEVATION_OFFSET.to_string(),
            "-resize",
            &(ELEVATION_TILE_SIZE * TILES_ACROSS_WIDTH / 4).to_string(),
            "+append",
        ];
        let top_row = &[
            "-size",
            &top_size,
            &format!("gray:{}", top_left.to_string_lossy()),
            &format!("gray:{}", top_right.to_string_lossy()),
        ];
        let bottom_row = &[
            "-size",
            &bottom_size,
            &format!("gray:{}", bottom_left.to_string_lossy()),
            &format!("gray:{}", bottom_right.to_string_lossy()),
        ];

        let crop_size = format!("{}x{}", ELEVATION_TILE_SIZE, ELEVATION_TILE_SIZE);
        let out_file = format!(
            "gray:{}",
            self.crops_dir().join("%[filename:tile]").to_string_lossy()
        );

        let output = Command::new("convert")
            .arg("-monitor")
            .args(&["-depth", "16"])
            .arg("(")
            .arg("-monitor")
            .args(top_row)
            .args(&offset_resize_and_append.clone())
            .arg(")")
            .arg("(")
            .arg("-monitor")
            .args(bottom_row)
            .args(&offset_resize_and_append.clone())
            .arg(")")
            .arg("-append")
            .args(&["-crop", &crop_size])
            .args(&["-set", "filename:tile", &name_template])
            .args(&["+repage", "+adjoin"])
            .arg(&out_file)
            .stderr(Stdio::inherit())
            .output()
            .chain_err(|| "Error while running `convert`")?;

        if !output.status.success() {
            return Err("`convert` returned with non-zero exit status".into());
        }

        Ok(())
    }

    fn create_crop_level(&self, level: u32) -> Result<()> {
        if self.crops_dir()
            .join(format!("{}_0_0.jpg", level))
            .is_file()
        {
            return Ok(());
        }

        self.pretty_log(
            3,
            level,
            MAX_LEVEL,
            &format!("Generating crops for level {}", level),
        );

        let num_crops_across_width = TILES_ACROSS_WIDTH / 2u32.pow(level);
        let num_crops_across_height = num_crops_across_width / 2;

        for crop_y in 0..num_crops_across_height {
            for crop_x in 0..num_crops_across_width {
                let left_x = crop_x * 2;
                let right_x = left_x + 1;
                let top_y = crop_y * 2;
                let bottom_y = top_y + 1;

                let top_left = format!("{}_{}_{}", level - 1, left_x, top_y);
                let top_right = format!("{}_{}_{}", level - 1, right_x, top_y);
                let bottom_left = format!("{}_{}_{}", level - 1, left_x, bottom_y);
                let bottom_right = format!("{}_{}_{}", level - 1, right_x, bottom_y);

                let elevation_crop_size =
                    format!("{}x{}", ELEVATION_TILE_SIZE, ELEVATION_TILE_SIZE);

                let resize_and_append = vec!["-resize", "50%", "+append"];

                let out_crop = format!("{}_{}_{}", level, crop_x, crop_y);
                let output = Command::new("convert")
                    .arg("(")
                    .arg(self.crops_dir().join(format!("{}.jpg", top_left)))
                    .arg(self.crops_dir().join(format!("{}.jpg", top_right)))
                    .args(resize_and_append.clone())
                    .arg(")")
                    .arg("(")
                    .arg(self.crops_dir().join(format!("{}.jpg", bottom_left)))
                    .arg(self.crops_dir().join(format!("{}.jpg", bottom_right)))
                    .args(resize_and_append.clone())
                    .arg(")")
                    .arg("-append")
                    .arg(self.crops_dir().join(format!("{}.jpg", out_crop)))
                    .stderr(Stdio::inherit())
                    .output()
                    .chain_err(|| "Error while running `convert`")?;

                if !output.status.success() {
                    return Err("`convert` returned with non-zero exit status".into());
                }

                let top_left = self.crops_dir().join(top_left).with_extension("elevation");
                let top_right = self.crops_dir().join(top_right).with_extension("elevation");
                let bottom_left = self.crops_dir().join(bottom_left).with_extension(
                    "elevation",
                );
                let bottom_right = self.crops_dir().join(bottom_right).with_extension(
                    "elevation",
                );
                let out_crop = self.crops_dir().join(out_crop).with_extension("elevation");

                let output = Command::new("convert")
                    .args(&["-depth", "16"])
                    .args(&["-size", &elevation_crop_size])
                    .arg("(")
                    .arg(format!("gray:{}", top_left.to_string_lossy()))
                    .arg(format!("gray:{}", top_right.to_string_lossy()))
                    .args(resize_and_append.clone())
                    .arg(")")
                    .arg("(")
                    .arg(format!("gray:{}", bottom_left.to_string_lossy()))
                    .arg(format!("gray:{}", bottom_right.to_string_lossy()))
                    .args(resize_and_append.clone())
                    .arg(")")
                    .arg("-append")
                    .arg(format!("gray:{}", out_crop.to_string_lossy()))
                    .stderr(Stdio::inherit())
                    .output()
                    .chain_err(|| "Error while running `convert`")?;

                if !output.status.success() {
                    return Err("`convert` returned with non-zero exit status".into());
                }
            }
        }

        Ok(())
    }

    fn create_tile_level(&self, level: u32) -> Result<()> {
        if self.tiles_dir()
            .join(format!("{}_0_0.elevation", level))
            .is_file()
        {
            return Ok(());
        }

        self.pretty_log(
            4,
            level + 1,
            MAX_LEVEL + 1,
            &format!("Generating tiles for level {}", level),
        );

        let num_crops_across_width = TILES_ACROSS_WIDTH / 2u32.pow(level);
        let num_crops_across_height = num_crops_across_width / 2;

        for tile_y in 0..num_crops_across_height {
            for tile_x in 0..num_crops_across_width {
                // For image tiles, no overlap is needed; just copy the relevant crop
                let nasa_file_name = format!("{}_{}_{}.jpg", level, tile_x, tile_y);
                std::fs::copy(
                    self.crops_dir().join(nasa_file_name.clone()),
                    self.tiles_dir().join(nasa_file_name.clone()),
                ).chain_err(|| "Error while copying Blue Marble crop to tiles directory")?;

                // For elevation tiles, we must crop and append tiles together
                let right_x = (tile_x + 1) % num_crops_across_width;
                let (below_y, gravity_below) = if tile_y == num_crops_across_height - 1 {
                    // For the last row, there's a special case. Instead of taking the top of the
                    // next row as usual, we use the bottom of the current row.
                    //
                    // In effect, this duplicates the bottom row of pixels of the bottom row of
                    // crops.
                    (tile_y, "southwest")
                } else {
                    (tile_y + 1, "northwest")
                };

                let elevation_crop_size =
                    format!("{}x{}", ELEVATION_TILE_SIZE, ELEVATION_TILE_SIZE);

                let top_left = format!("{}_{}_{}.elevation", level, tile_x, tile_y);
                let top_right = format!("{}_{}_{}.elevation", level, right_x, tile_y);
                let bottom_left = format!("{}_{}_{}.elevation", level, tile_x, below_y);
                let bottom_right = format!("{}_{}_{}.elevation", level, right_x, below_y);

                let top_left =
                    format!("gray:{}", self.crops_dir().join(top_left).to_string_lossy());
                let top_right = format!(
                    "gray:{}",
                    self.crops_dir().join(top_right).to_string_lossy()
                );
                let bottom_left = format!(
                    "gray:{}",
                    self.crops_dir().join(bottom_left).to_string_lossy()
                );
                let bottom_right = format!(
                    "gray:{}",
                    self.crops_dir().join(bottom_right).to_string_lossy()
                );

                let top_right_crop = format!("1x{}+0+0", ELEVATION_TILE_SIZE);
                let bottom_left_crop = format!("{}x1+0+0", ELEVATION_TILE_SIZE);
                let bottom_right_crop = format!("1x1+0+0");

                let out_tile = format!("{}_{}_{}.elevation", level, tile_x, tile_y);
                let out_tile =
                    format!("gray:{}", self.tiles_dir().join(out_tile).to_string_lossy());

                let top_left = ["(", &top_left, ")"];
                let top_right = ["(", &top_right, "-crop", &top_right_crop, ")"];
                let bottom_left = [
                    "(",
                    &bottom_left,
                    "-gravity",
                    gravity_below,
                    "-crop",
                    &bottom_left_crop,
                    ")",
                ];
                let bottom_right = [
                    "(",
                    &bottom_right,
                    "-gravity",
                    gravity_below,
                    "-crop",
                    &bottom_right_crop,
                    ")",
                ];

                let output = Command::new("convert")
                    .args(&["-depth", "16"])
                    .args(&["-size", &elevation_crop_size])
                    .arg("(")
                    .args(&top_left)
                    .args(&top_right)
                    .arg("+append")
                    .arg(")")
                    .arg("(")
                    .args(&bottom_left)
                    .args(&bottom_right)
                    .arg("+append")
                    .arg(")")
                    .arg("-append")
                    .args(&["-format", "%[fx:minima] %[fx:maxima]"])
                    .args(&["-write", "info:-"])
                    .arg(out_tile)
                    .stderr(Stdio::inherit())
                    .output()
                    .chain_err(|| "Error while running `convert`")?;

                if !output.status.success() {
                    return Err("`convert` returned with non-zero exit status".into());
                }

                let output = String::from_utf8(output.stdout).unwrap();
                let output_parts: Vec<_> = output.split(" ").collect();
                let min = f32::from_str(output_parts[0]).unwrap();
                let max = f32::from_str(output_parts[1]).unwrap();

                let min = (min * u16::max_value() as f32) as u16;
                let max = (max * u16::max_value() as f32) as u16;

                let metadata_path = self.tiles_dir().join(format!(
                    "{}_{}_{}.json",
                    level,
                    tile_x,
                    tile_y
                ));
                let mut metadata_file = File::create(metadata_path).chain_err(
                    || "Unable to create first pass tile metadata file",
                )?;

                let metadata = FirstPassTileMetadata {
                    min_elevation: min.saturating_sub(ELEVATION_OFFSET),
                    max_elevation: max.saturating_sub(ELEVATION_OFFSET),
                };

                let to_write = serde_json::to_string(&metadata).unwrap();
                metadata_file.write_all(to_write.as_bytes()).chain_err(
                    || "Error writing out first pass tile metadata",
                )?;
            }
        }

        Ok(())
    }

    fn create_polygon_data(&self) -> Result<()> {
        // if self.output_dir.join("polygons.json").is_file() {
        //     return Ok(())
        // }

        let mut features_file = File::open(&self.features_file).chain_err(
            || "Could not open features file",
        )?;
        let mut geojson = String::new();
        features_file.read_to_string(&mut geojson).chain_err(
            || "Could not read from features file",
        )?;

        let geojson: GeoJson = geojson.parse().unwrap();
        let feature_collection = match geojson {
            GeoJson::FeatureCollection(fc) => fc,
            _ => panic!("Unexpected geojson object type!"),
        };

        let mut polygons = Vec::new();
        let mut polygon_properties = Vec::new();

        for feature in feature_collection.features {
            let properties = feature.properties.unwrap_or(serde_json::Map::new());
            let geometry: geo::Geometry<f32> = feature.geometry.unwrap().value.try_into().unwrap();

            match geometry {
                geo::Geometry::MultiPolygon(mp) => {
                    for polygon in mp {
                        polygons.push(polygon);
                        polygon_properties.push(properties.clone());
                    }
                }
                geo::Geometry::Polygon(p) => {
                    polygons.push(p);
                    polygon_properties.push(properties);
                }
                _ => panic!(),
            }
        }

        self.create_final_tile_metadatas(&polygons)?;
        self.create_polygon_file(&polygons, &polygon_properties)?;
        Ok(())
    }

    fn create_final_tile_metadatas(&self, polygons: &[geo::Polygon<f32>]) -> Result<()> {
        let mut tile_metadatas: BTreeMap<(u32, u32, u32), Vec<u64>> = BTreeMap::new();
        self.pretty_log(5, 1, 1, "Generating tile metadata files");

        for level in 0..MAX_LEVEL + 1 {
            let num_tiles_across_width = TILES_ACROSS_WIDTH / 2u32.pow(level);
            let num_tiles_across_height = num_tiles_across_width / 2;

            for (polygon_index, polygon) in polygons.iter().enumerate() {
                let polygon_bbox = polygon.exterior.bbox().unwrap();
                let x_min = (polygon_bbox.xmin + 180.0) / 360.0;
                let x_max = (polygon_bbox.xmax + 180.0) / 360.0;
                let y_min = (-polygon_bbox.ymin + 90.0) / 180.0;
                let y_max = (-polygon_bbox.ymax + 90.0) / 180.0;

                let xs = (x_min * num_tiles_across_width as f32) as u32..
                    (x_max * num_tiles_across_width as f32) as u32 + 1;
                let ys = (y_max * num_tiles_across_height as f32) as u32..
                    (y_min * num_tiles_across_height as f32) as u32 + 1;

                for x in xs {
                    for y in ys.clone() {
                        tile_metadatas
                            .entry((level, x, y))
                            .or_insert(Vec::new())
                            .push(polygon_index as u64);
                    }
                }
            }
        }

        for level in 0..MAX_LEVEL + 1 {
            let num_tiles_across_width = TILES_ACROSS_WIDTH / 2u32.pow(level);
            let num_tiles_across_height = num_tiles_across_width / 2;

            for x in 0..num_tiles_across_width {
                for y in 0..num_tiles_across_height {
                    let metadata_path =
                        self.tiles_dir().join(format!("{}_{}_{}.json", level, x, y));
                    let metadata_file = File::open(&metadata_path).chain_err(
                        || "Unable to read in first pass metadata",
                    )?;
                    let first_pass_metadata: FirstPassTileMetadata =
                        serde_json::from_reader(metadata_file).chain_err(
                            || "Error parsing first pass metadata",
                        )?;

                    let default_polygons = Vec::new();
                    let polygons = tile_metadatas.get(&(level, x, y)).unwrap_or(
                        &default_polygons,
                    );

                    let tile_metadata = TileMetadata {
                        min_elevation: first_pass_metadata.min_elevation,
                        max_elevation: first_pass_metadata.max_elevation,
                        polygons: polygons.clone(),
                    };

                    let mut metadata_file = File::create(metadata_path).chain_err(
                        || "Unable to open metadata file",
                    )?;

                    let to_write = serde_json::to_string(&tile_metadata).unwrap();
                    metadata_file.write_all(to_write.as_bytes()).chain_err(
                        || "Error writing out tile metadata",
                    )?;
                }
            }
        }

        Ok(())
    }

    fn create_polygon_file(
        &self,
        polygons: &[geo::Polygon<f32>],
        polygon_properties: &[serde_json::Map<String, serde_json::Value>],
    ) -> Result<()> {
        self.pretty_log(6, 1, 1, "Generating polygon point data file");

        let multi_level_polygons = polygons
            .iter()
            .zip(polygon_properties)
            .map(|(polygon, properties)| {
                let levels = (0..MAX_LEVEL + 1)
                    .map(|level| {
                        let polygon =
                            polygon.simplify(&self.simplification_epsilons[level as usize]);
                        let mut points = Vec::new();

                        points.extend(polygon.exterior.into_iter().map(
                            |point| (point.x(), point.y()),
                        ));
                        for ring in polygon.interiors {
                            points.extend(ring.into_iter().map(|point| (point.x(), point.y())));
                        }

                        points
                    })
                    .collect();

                MultiLevelPolygon {
                    properties: properties.clone(),
                    levels: levels,
                }
            })
            .collect();

        let polygon_path = self.output_dir.join("polygons.json");
        let mut polygon_file = File::create(polygon_path).chain_err(
            || "Unable to open polygons file",
        )?;

        let polygon_point_data = PolygonPointData { polygons: multi_level_polygons };

        let to_write = serde_json::to_string(&polygon_point_data).unwrap();
        polygon_file.write_all(to_write.as_bytes()).chain_err(
            || "Error writing out polygon point data",
        )?;

        Ok(())
    }

    fn pretty_log(&self, step: u32, sub_step: u32, num_sub_steps: u32, message: &str) {
        let to_print = format!(
            "[{}/{}] [{}/{}]: {}",
            step,
            NUM_STEPS,
            sub_step,
            num_sub_steps,
            message
        );

        println!("{}", to_print.bold().green());
    }

    fn crops_dir(&self) -> PathBuf {
        self.output_dir.join("crops")
    }

    fn tiles_dir(&self) -> PathBuf {
        self.output_dir.join("tiles")
    }
}
