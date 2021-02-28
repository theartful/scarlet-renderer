#![allow(dead_code)]

use rand::distributions::{Distribution, Uniform};
use scarlet::camera::*;
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
    let image_width: u32 = 400;
    let image_height: u32 = (image_width as Float / aspect_ratio) as u32;

    let focal_length = 1.0;
    let fov = Float::PI * 0.26;

    let camera = PerspectiveCamera::init(
        Transformf::look_at(
            Point3::init(0., 0., 2.),
            Point3::init(0., 0., -1.),
            Vector3f::init(0., 1., 0.),
        ),
        focal_length,
        fov,
        focal_length / 3.,
    );

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

    let du = 1. / (image_height - 1) as Float;
    let dv = 1. / (image_width - 1) as Float;
    let uniform_du = Uniform::new(-du / 2., du / 2.);
    let uniform_dv = Uniform::new(-dv / 2., dv / 2.);
    let mut rng = rand::thread_rng();

    for jj in 0..image_height {
        let j = image_height - jj - 1;
        for i in 0..image_width {
            let mut color = Vector3f::new();
            for _ in 0..sample_per_pixel {
                let u = 2. * i as Float / (image_width as Float - 1.) - 1.
                    + uniform_du.sample(&mut rng);
                let v = 2. * j as Float / (aspect_ratio * (image_height as Float - 1.))
                    - 1. / aspect_ratio
                    + uniform_dv.sample(&mut rng);
                let ray =
                    camera.generate_ray(CameraSample::init(Point2f::init(u, v), Point2f::new()));
                color += ray_color(Rayfi::from(ray), &primitive_list, 0);
            }
            write_color(color / sample_per_pixel as Float);
        }
    }
}
