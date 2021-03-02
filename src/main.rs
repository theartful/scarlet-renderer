#![allow(dead_code)]

use rand::distributions::{Distribution, Uniform};
use scarlet::camera::*;
use scarlet::interpolation::*;
use scarlet::primitive::*;
use scarlet::ray::*;
use scarlet::sphere::*;
use scarlet::vector::*;
use std::ops::DerefMut;

use rand::SeedableRng;

static mut rng_box: Option<Box<rand::rngs::StdRng>> = None;

fn get_rng() -> &'static mut rand::rngs::StdRng {
    unsafe {
        match &mut rng_box {
            Some(rng) => rng.deref_mut(),
            None => panic!(),
        }
    }
}

fn init_rng() {
    unsafe {
        rng_box = Some(Box::new(rand::rngs::StdRng::from_seed([
            0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
            24, 25, 26, 27, 28, 29, 30, 31,
        ])));
    }
}

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
    let rng = get_rng();

    let theta: Float = uniform_theta.sample(rng);
    let phi: Float = uniform_phi.sample(rng);

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
            let n = surface_interaction.n;
            let p = surface_interaction.p;
            let u = surface_interaction.dpdu.normalize();
            let v = surface_interaction.dpdv.normalize();
            let dir = Vector3fi::from(sample_direction(n, u, v));

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
    init_rng();
    let rng = get_rng();

    let aspect_ratio: Float = 16. / 9.;
    let image_width: u32 = 1000;
    let image_height: u32 = (image_width as Float / aspect_ratio) as u32;

    let focus_distance = 1.;
    let lens_radius = 0.01;
    let fov = Float::PI * 0.26;

    let camera = ThinLensCamera::init(
        Transformf::look_at(
            Point3::init(0., 0., 2.),
            Point3::init(0., 0., -1.),
            Vector3f::init(0., 1., 0.),
        ),
        focus_distance,
        fov,
        lens_radius,
        0.01,
        1000.,
    );

    println!("P3");
    println!("{} {}", image_width, image_height);
    println!("255");

    let sphere = Sphere::init(0.5, -1., 1., 2. * Float::PI);
    let transform = Transform::translate(Vector3f::init(0., 0., -1.));
    let primitive = TransformedPrimitive::init(sphere, transform);

    let sphere2 = Sphere::init(100., -200., 200., 2. * Float::PI);
    let transform2 = Transform::translate(Vector3f::init(0., -100.5, -1.));
    let primitive2 = TransformedPrimitive::init(sphere2, transform2);

    let mut primitive_list = PrimitiveList::new();
    primitive_list.add_primitive(primitive2);
    primitive_list.add_primitive(primitive);

    let sample_per_pixel = 100;

    let du = 1. / (image_height - 1) as Float;
    let dv = 1. / (image_width - 1) as Float;
    let uniform_du = Uniform::new(0., du);
    let uniform_dv = Uniform::new(0., dv);
    let uniform_aperture = Uniform::new(-1., 1.);

    for jj in 0..image_height {
        let j = image_height - jj - 1;
        for i in 0..image_width {
            let mut color = Vector3f::new();
            for _ in 0..sample_per_pixel {
                let u = 2. * i as Float / (image_width as Float - 1.) - 1. + uniform_du.sample(rng);
                let v = 2. * j as Float / (aspect_ratio * (image_height as Float - 1.))
                    - 1. / aspect_ratio
                    + uniform_dv.sample(rng);
                let ray = camera.generate_ray(CameraSample::init(
                    Point2f::init(u, v),
                    Point2f::init(uniform_aperture.sample(rng), uniform_aperture.sample(rng)),
                ));
                color += ray_color(Rayfi::from(ray), &primitive_list, 0);
            }
            write_color(color / sample_per_pixel as Float);
        }
    }

    // let i = 499;
    // let j = 281;
    // let u = 2. * i as Float / (image_width as Float - 1.) - 1.;
    // let v = 2. * j as Float / (aspect_ratio * (image_height as Float - 1.)) - 1. / aspect_ratio;
    // let mut ray =
    //     Rayfi::from(camera.generate_ray(CameraSample::init(Point2f::init(u, v), Point2f::new())));
    // match primitive_list.intersecti(&mut ray) {
    //     Some(_) => {
    //         eprintln!("OK OK!\n");
    //     }
    //     None => {
    //         eprintln!("NOT OK AT ALL!\n");
    //     }
    // }
}
