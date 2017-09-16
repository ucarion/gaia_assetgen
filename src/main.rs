extern crate gaia_assetgen;

use gaia_assetgen::PrepareAssetsTask;

fn main() {
    PrepareAssetsTask::new()
        .with_noaa_globe_dir("assets/noaa_globe".into())
        .with_nasa_blue_marble_dir("assets/nasa_blue_marble".into())
        .with_features_file("assets/ne_10m_admin_0_countries.geojson".into())
        .with_output_dir("assets/generated".into())
        .run()
        .unwrap();
}
