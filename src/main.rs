use macroquad::prelude::*;
use std::any::Any;

struct Sphere {
    radius: f32,
    color: Color
}

enum Shape {
    Sphere
}

trait Shaped {
    fn shape_type(&self) -> Shape;
    fn get_centre_of_mass(&self) -> Vec3;
    fn radius(&self) -> f32;
    fn as_any(&self) -> &dyn Any;
    fn draw(&self, pos:Vec3);
}

impl Shaped for Sphere {
    fn shape_type(&self) -> Shape {
        Shape::Sphere
    }

    fn get_centre_of_mass(&self) -> Vec3 {
        Vec3::ZERO
    }

    fn radius(&self) -> f32 {
        self.radius
    }

    fn as_any(&self) -> &dyn Any {
        self
    }

    fn draw(&self, pos: Vec3) {
        draw_sphere(pos, self.radius, None, self.color);
    }
}

struct Body {
    position: Vec3,
    orientation: Quat,
    linear_veclocity: Vec3,
    inv_mass: f32,
    elasticity: f32,
    shape: Box<dyn Shaped>
}

impl Body {
    fn get_centre_of_mass_world_space(&self) -> Vec3 {
        let com = self.shape.get_centre_of_mass();
        self.position + self.orientation.mul_vec3(com) 
    }
    fn get_centre_of_mass_model_space(&self) -> Vec3 {
        self.shape.get_centre_of_mass()
    }
    fn world_space_to_body_space(&self, world_pt:Vec3) -> Vec3 {
        self.orientation.inverse().mul_vec3(
            world_pt - self.get_centre_of_mass_world_space())
    }
    fn body_space_to_world_space(&self, body_pt:Vec3) -> Vec3 {
        self.get_centre_of_mass_world_space() +
            self.orientation.mul_vec3(body_pt)
    }
    fn apply_impluse_linear(&mut self, impulse:Vec3) {
        if self.inv_mass == 0. {
            return;
        }
        self.linear_veclocity += impulse * self.inv_mass;
    }
    fn radius(&self) -> f32 {
        self.shape.radius()
    }
    fn draw(&self) {
        self.shape.draw(self.position);
    }
}

struct Scene {
    bodies: Vec<Body>
}

struct Contact {
    pt_on_a_world_space: Vec3,
    pt_on_b_world_space: Vec3,
    pt_on_a_local_space: Vec3,
    pt_on_b_local_space: Vec3,
    normal: Vec3, 
    separation_distance: f32,
    time_of_impact: f32,
    body_a: usize,
    body_b: usize
}

fn intersect(bodies: &[Body], i: usize, j: usize) -> Option<Contact> {
    let (a, b) = (&bodies[i], &bodies[j]);
    let ab = b.position - a.position;
    let radius_ab = b.radius() + a.radius();
    if ab.length_squared() <= radius_ab * radius_ab {
        let normal = ab.normalize();
        let pt_on_a_world_space = a.position + normal * a.radius();
        let pt_on_b_world_space = b.position - normal * b.radius();
        Some(Contact {
            pt_on_a_world_space,
            pt_on_b_world_space,
            pt_on_a_local_space: Vec3::ZERO,
            pt_on_b_local_space: Vec3::ZERO,
            normal, 
            separation_distance: 0.,
            time_of_impact: 0.,
            body_a: i,
            body_b: j
        })
    } else {
        None
    }
}

fn resolve_contact(bodies: &mut[Body], contact: &Contact) {
    let vec_impulse_j:Vec3;
    {
        let (a, b) = (&bodies[contact.body_a], &bodies[contact.body_b]);
        let vab = a.linear_veclocity - b.linear_veclocity;
        let elasticity = a.elasticity * b.elasticity;
        let impulse_j = -(1. + elasticity) * vab.dot(contact.normal) / (a.inv_mass + b.inv_mass);
        vec_impulse_j = contact.normal * impulse_j;
    }
    bodies[contact.body_a].apply_impluse_linear(vec_impulse_j);
    bodies[contact.body_b].apply_impluse_linear(vec_impulse_j * -1.);

    // move colliders to just outside each other but keeping combined centre of
    // mass constant
    let (a, b) = (&bodies[contact.body_a], &bodies[contact.body_b]);
    let ta = a.inv_mass / (a.inv_mass + b.inv_mass);
    let tb = b.inv_mass / (a.inv_mass + b.inv_mass);
    let ds = contact.pt_on_b_world_space - contact.pt_on_a_world_space;
    bodies[contact.body_a].position += ds * ta;
    bodies[contact.body_b].position -= ds * tb;
}

impl Scene {
    fn new() -> Scene {
        let mut bodies = Vec::<Body>::new();
        bodies.push(
            Body {
                position: Vec3::new(0., 10., 0.),
                orientation: Quat::IDENTITY,
                linear_veclocity: Vec3::ZERO,
                inv_mass: 1.,
                elasticity: 0.5,
                shape: Box::new(Sphere{radius:1., color:BLUE})
            }
        );
        bodies.push(
            Body {
                position: Vec3::new(0., -1000., 0.),
                orientation: Quat::IDENTITY,
                linear_veclocity: Vec3::ZERO,
                inv_mass: 0.,
                elasticity: 1.,
                shape: Box::new(Sphere{radius:1000., color:GREEN})
            }
        );        Scene {
            bodies
        }
    }

    fn update(&mut self, dt_sec: f32) {
        for body in &mut self.bodies {
            let mass = 1. / body.inv_mass;
            let impulse_gravity = Vec3::new(0., -10., 0.) * mass * dt_sec;
            body.apply_impluse_linear(impulse_gravity);
        }
        let n = self.bodies.len();
        for i in 0..n {
            for j in i + 1..n {
                let (a, b) = (&self.bodies[i], &self.bodies[j]);
                if a.inv_mass == 0. && b.inv_mass == 0. {
                    continue
                }
                if let Some(contact) = intersect(&self.bodies[..], i, j) {
                    resolve_contact(&mut self.bodies[..], &contact);
                }
            }
        }

        for body in &mut self.bodies {
            body.position += body.linear_veclocity * dt_sec;
        }
   }

    fn draw(&self) {
        for body in &self.bodies {
            body.draw();
        }
    }
}

#[macroquad::main("3D")]
async fn main() {
    let mut scene = Scene::new();

    loop {
        scene.update(1.0 / 60.0);
        clear_background(LIGHTGRAY);

        // Going 3d!

        set_camera(&Camera3D {
            position: vec3(-20., 15., 0.),
            up: vec3(0., 1., 0.),
            target: vec3(0., 0., 0.),
            ..Default::default()
        });

        draw_grid(20, 1., BLACK, GRAY);

        scene.draw();
        // Back to screen space, render some text

        set_default_camera();
        draw_text("WELCOME TO 3D WORLD", 10.0, 20.0, 30.0, BLACK);

        next_frame().await
    }
}