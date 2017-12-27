use std::path::Path;
use std::process::{Command, Stdio};
use std::str::FromStr;

use errors::*;

#[derive(Debug)]
pub struct Convert {
    command: Command,
}

impl Convert {
    pub fn new() -> Convert {
        let mut command = Command::new("convert");
        command.stderr(Stdio::inherit());

        Convert { command }
    }

    pub fn monitor(&mut self) -> &mut Convert {
        self.command.arg("-monitor");
        self
    }

    pub fn grayscale_input(&mut self, size: (u32, u32), path: &Path) -> &mut Convert {
        let size_arg = format!("{}x{}", size.0, size.1);
        self.command.args(&["-depth", "16", "-size", &size_arg]);
        self.command.arg(format!("gray:{}", path.to_string_lossy()));
        self
    }

    pub fn input(&mut self, path: &Path) -> &mut Convert {
        self.command.arg(path);
        self
    }

    pub fn grayscale_output(&mut self, path: &Path) -> &mut Convert {
        self.command.arg("+adjoin");
        self.command.args(&["-depth", "16"]);
        self.command.arg(format!("gray:{}", path.to_string_lossy()));
        self
    }

    pub fn output(&mut self, path: &Path) -> &mut Convert {
        self.command.arg(path);
        self
    }

    pub fn group<F>(&mut self, f: F) -> &mut Convert
    where
        F: FnOnce(&mut Convert) -> &mut Convert,
    {
        self.command.arg("(");
        f(self);
        self.command.arg(")");
        self
    }

    pub fn crops(&mut self, size: u32) -> &mut Convert {
        let crop_fmt = format!("{}x{}", size, size);

        self.command.args(&["-crop", &crop_fmt]);
        self
    }

    pub fn crop_one(&mut self, size: (u32, u32), offset: (u32, u32)) -> &mut Convert {
        let crop_fmt = format!("{}x{}+{}+{}", size.0, size.1, offset.0, offset.1);

        self.command.args(&["-crop", &crop_fmt]);
        self
    }

    pub fn append_horizontally(&mut self) -> &mut Convert {
        self.command.arg("+append");
        self
    }

    pub fn append_vertically(&mut self) -> &mut Convert {
        self.command.arg("-append");
        self
    }

    pub fn resize(&mut self, resize: &str) -> &mut Convert {
        self.command.args(&["-resize", resize]);
        self
    }

    pub fn offset_each_pixel(&mut self, offset: u16) -> &mut Convert {
        self.command.args(
            &[
                "-evaluate",
                "addmodulus",
                &offset.to_string(),
            ],
        );
        self
    }

    pub fn report_max_min(&mut self) -> &mut Convert {
        self.command.args(
            &["-format", "%[fx:minima] %[fx:maxima]", "-write", "info:-"],
        );
        self
    }

    pub fn run(&mut self) -> Result<String> {
        let output = self.command.output().chain_err(
            || "Error when running `convert`",
        )?;

        if !output.status.success() {
            return Err("`convert` returned with non-zero exit status".into());
        }

        String::from_utf8(output.stdout).chain_err(|| "Error parsing output as UTF-8")
    }

    pub fn run_with_max_min(&mut self) -> Result<(f32, f32)> {
        let output = self.run()?;

        let output_parts: Vec<_> = output.split(" ").collect();
        let min = f32::from_str(output_parts[0]).unwrap();
        let max = f32::from_str(output_parts[1]).unwrap();
        Ok((min, max))
    }
}