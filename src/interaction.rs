use crate::transform::*;
use crate::vector::*;

pub struct SurfaceInteraction {
    pub p: Point3fi,
    pub n: Normal3f,
    pub uv: Point2f,
    pub dpdu: Vector3f,
    pub dpdv: Vector3f,
    // dndu: Normal3f,
    // dndv: Normal3f,
}

impl SurfaceInteraction {
    pub fn init(p: Point3fi, uv: Point2f, dpdu: Vector3f, dpdv: Vector3f) -> Self {
        SurfaceInteraction {
            p,
            n: Normal3f::from(dpdu.cross(dpdv).normalize()),
            uv,
            dpdu,
            dpdv,
        }
    }
}

impl TransformMap<SurfaceInteraction> for Transformf {
    type Output = SurfaceInteraction;

    fn map(&self, interaction: SurfaceInteraction) -> Self::Output {
        SurfaceInteraction::init(
            self.map(interaction.p),
            interaction.uv,
            self.map(interaction.dpdu),
            self.map(interaction.dpdv),
        )
    }
}
