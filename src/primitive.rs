pub use crate::bbox::*;
pub use crate::interaction::*;
pub use crate::ray::*;
pub use crate::transform::*;
pub use crate::vector::*;
pub use crate::vector::*;
use std::rc::Rc;

pub trait Primitive: Bounded3<Float> {
    fn area(&self) -> Float;

    // one of intersect/intersecti has to be implemented
    fn intersect(&self, ray: &mut Rayf) -> Option<SurfaceInteraction> {
        let mut rayfi = Rayfi::from(*ray);
        let res = self.intersecti(&mut rayfi);
        ray.tmax = rayfi.tmax.approx();
        res
    }
    fn intersecti(&self, ray: &mut Rayfi) -> Option<SurfaceInteraction> {
        let mut rayf = Rayf::from(*ray);
        let res = self.intersect(&mut rayf);
        ray.tmax = Intervalf::from(rayf.tmax);
        res
    }
    fn fast_intersect(&self, ray: Rayf) -> bool {
        self.fast_intersecti(Rayfi::from(ray))
    }
    fn fast_intersecti(&self, ray: Rayfi) -> bool {
        self.fast_intersect(Rayf::from(ray))
    }
}

pub type PrimitiveHandle = Rc<dyn Primitive>;
pub fn alloc_primitive<T: Primitive + 'static>(primitive: T) -> PrimitiveHandle {
    Rc::new(primitive) as PrimitiveHandle
}

pub struct TransformedPrimitive {
    primitive_handle: Rc<dyn Primitive>,
    otow: Transformf,
}

impl TransformedPrimitive {
    pub fn init<T: Primitive + 'static>(primitive: T, otow: Transformf) -> Self {
        TransformedPrimitive {
            primitive_handle: alloc_primitive(primitive),
            otow,
        }
    }
    fn object_to_world<T, U>(&self, t: T) -> U
    where
        Transformf: TransformMap<T, Output = U>,
    {
        self.otow.map(t)
    }
    fn world_to_object<T, U>(&self, t: T) -> U
    where
        Transformf: TransformMap<T, Output = U>,
    {
        self.otow.map_inverse(t)
    }
}

impl Bounded3<Float> for TransformedPrimitive {
    fn bbox(&self) -> Bbox3f {
        self.object_to_world((*self.primitive_handle).bbox())
    }
}

impl Primitive for TransformedPrimitive {
    fn intersecti(&self, ray: &mut Rayfi) -> Option<SurfaceInteraction> {
        // note that "t" computed in object space is the same as in world space
        let mut transformed_ray = self.world_to_object(*ray);
        match (*self.primitive_handle).intersecti(&mut transformed_ray) {
            None => None,
            Some(surface_interaction) => {
                ray.tmax = transformed_ray.tmax;
                Some(self.object_to_world(surface_interaction))
            }
        }
    }
    fn fast_intersecti(&self, ray: Rayfi) -> bool {
        (*self.primitive_handle).fast_intersecti(self.world_to_object(ray))
    }
    fn area(&self) -> Float {
        // TODO: how is the area mapped under linear transformation with sheer and scaling effects?
        (*self.primitive_handle).area()
    }
}

impl Bounded3<Float> for PrimitiveHandle {
    fn bbox(&self) -> Bbox3<Float> {
        (**self).bbox()
    }
}

pub struct PrimitiveList {
    primitives: Vec<PrimitiveHandle>,
}

impl PrimitiveList {
    pub fn new() -> PrimitiveList {
        PrimitiveList {
            primitives: Vec::new(),
        }
    }

    pub fn add_primitive<T: Primitive + 'static>(&mut self, primitive: T) {
        self.primitives.push(alloc_primitive(primitive));
    }
}

impl Bounded3<Float> for PrimitiveList {
    fn bbox(&self) -> Bbox3f {
        let mut bbox = Bbox3f::new();
        for primitive in self.primitives.iter() {
            bbox += (*primitive).bbox();
        }
        bbox
    }
}

impl Primitive for PrimitiveList {
    fn intersect(&self, ray: &mut Rayf) -> Option<SurfaceInteraction> {
        let mut result: Option<SurfaceInteraction> = None;
        for primitive_handle in self.primitives.iter() {
            match (*primitive_handle).intersect(ray) {
                Some(interaction) => result = Some(interaction),
                None => {}
            }
        }
        result
    }
    fn intersecti(&self, ray: &mut Rayfi) -> Option<SurfaceInteraction> {
        let mut result: Option<SurfaceInteraction> = None;
        for primitive_handle in self.primitives.iter() {
            match (*primitive_handle).intersecti(ray) {
                Some(interaction) => result = Some(interaction),
                None => {}
            }
        }
        result
    }
    fn fast_intersect(&self, ray: Rayf) -> bool {
        for primitive_handle in self.primitives.iter() {
            if (*primitive_handle).fast_intersect(ray) {
                return true;
            }
        }
        return false;
    }
    fn fast_intersecti(&self, ray: Rayfi) -> bool {
        for primitive_handle in self.primitives.iter() {
            if (*primitive_handle).fast_intersecti(ray) {
                return true;
            }
        }
        return false;
    }
    fn area(&self) -> Float {
        let mut tot_area: Float = Float::zero();
        for primitive_handle in self.primitives.iter() {
            tot_area += (*primitive_handle).area();
        }
        tot_area
    }
}
