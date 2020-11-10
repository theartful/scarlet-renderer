#![allow(dead_code)]

use rand::distributions::{Distribution, Uniform};
use scarlet::interpolation::*;
use scarlet::primitive::*;
use scarlet::ray::*;
use scarlet::sphere::*;
use scarlet::vector::*;

fn write_color(color: Vector3f) {
    print!(
        "{} {} {}\n",
        (255. * color.x).round(),
        (255. * color.y).round(),
        (255. * color.z).round()
    )
}

fn sample_direction(n: Normal3f, u: Vector3f, v: Vector3f) -> Vector3f {
    let uniform_theta = Uniform::new(0., Float::PI / 2.);
    let uniform_phi = Uniform::new(0., Float::PI * 2.);
    let mut rng = rand::thread_rng();

    let theta: Float = uniform_theta.sample(&mut rng);
    let phi: Float = uniform_phi.sample(&mut rng);

    let (sin_theta, cos_theta) = theta.sin_cos();
    let (sin_phi, cos_phi) = phi.sin_cos();

    (Vector3f::from(n) * (1. + cos_theta) + u * sin_theta * cos_phi + v * sin_theta * sin_phi)
        .normalize()
}

fn ray_color(ray: Rayfi, primitive: &impl Primitive, depth: u32) -> Vector3f {
    let max_depth = 50;
    if depth >= max_depth {
        return Vector3f::new();
        // eprintln!("depth = {}", depth);
    }

    let mut mut_ray = ray;

    match primitive.intersecti(&mut mut_ray) {
        Some(surface_interaction) => {
            let n = Vector3f::from(surface_interaction.n);
            let p = surface_interaction.p;
            let u = surface_interaction.dpdu.normalize();
            let v = surface_interaction.dpdv.normalize();
            let dir = Vector3fi::from(sample_direction(Normal3f::from(n), u, v));

            ray_color(
                Rayfi::init(p, dir, Intervalf::highest()),
                primitive,
                depth + 1,
            ) * 0.5
        }
        None => lerp(
            Vector3f::init(1., 1., 1.),
            Vector3f::init(0.5, 0.7, 1.),
            (mut_ray.direction.y.approx() + 1.0) / 2.0,
        ),
    }
}

fn main() {
    let aspect_ratio: Float = 16. / 9.;
    let image_width: u32 = 1000;
    let image_height: u32 = (image_width as Float / aspect_ratio) as u32;

    let viewport_height: Float = 2.;
    let viewport_width: Float = aspect_ratio * viewport_height;
    let focal_length = 1.0;

    let origin = Point3f::init(0., 0., 0.);
    let horizontal = Vector3f::init(viewport_width, 0., 0.);
    let vertical = Vector3f::init(0., viewport_height, 0.);
    let lower_left_corner =
        origin - horizontal / 2. - vertical / 2. - Vector3f::init(0., 0., focal_length);

    println!("P3");
    println!("{} {}", image_width, image_height);
    println!("255");

    let sphere = Sphere::init(0.5, -0.25, 0.25, 2. * Float::PI);
    let transform = Transform::translate(Vector3f::init(0., 0., -1.));
    let primitive = TransformedPrimitive::init(sphere, transform);

    let sphere2 = Sphere::init(100., -100., 100., 2. * Float::PI);
    let transform2 = Transform::translate(Vector3f::init(0., -100.5, -1.));
    let primitive2 = TransformedPrimitive::init(sphere2, transform2);

    let mut primitive_list = PrimitiveList::new();
    primitive_list.add_primitive(primitive2);
    primitive_list.add_primitive(primitive);

    let sample_per_pixel = 100;

    for jj in 0..image_height {
        let j = image_height - jj - 1;
        for i in 0..image_width {
            let du = 1. / (image_height - 1) as Float;
            let dv = 1. / (image_width - 1) as Float;
            let mut color = Vector3f::new();
            for s in 0..sample_per_pixel {
                let v = (j as Float / (image_height - 1) as Float)
                    + (s as Float / (sample_per_pixel - 1) as Float) * du;
                let u = (i as Float / (image_width - 1) as Float)
                    + (s as Float / (sample_per_pixel - 1) as Float) * dv;
                let dir = (lower_left_corner + horizontal * u + vertical * v - origin).normalize();
                let r = Rayf::init(origin, dir, Float::highest());
                color += ray_color(Rayfi::from(r), &primitive_list, 0);
            }
            write_color(color / sample_per_pixel as Float);
        }
    }
}
