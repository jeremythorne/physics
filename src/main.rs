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
    fn inertial_tensor(&self) -> Mat3;
    fn radius(&self) -> f32;
    fn as_any(&self) -> &dyn Any;
    fn draw(&self, pos:Vec3, texture:Texture2D);
}

impl Shaped for Sphere {
    fn shape_type(&self) -> Shape {
        Shape::Sphere
    }

    fn get_centre_of_mass(&self) -> Vec3 {
        Vec3::ZERO
    }

    fn inertial_tensor(&self) -> Mat3 {
        let f = 2. * self.radius * self.radius() / 5.;
        Mat3::from_diagonal(Vec3::new(f, f, f))
    }

    fn radius(&self) -> f32 {
        self.radius
    }

    fn as_any(&self) -> &dyn Any {
        self
    }

    fn draw(&self, pos: Vec3, texture:Texture2D) {
        draw_sphere(pos, self.radius, texture, self.color);
    }
}

struct Body {
    position: Vec3,
    orientation: Quat,
    linear_veclocity: Vec3,
    angular_veclocity: Vec3,
    inv_mass: f32,
    elasticity: f32,
    friction: f32,
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
    fn get_inverse_inertial_tensor_body_space(&self) -> Mat3 {
        self.shape.inertial_tensor().inverse() * self.inv_mass
    }
    fn get_inverse_inertial_tensor_world_space(&self) -> Mat3 {
        let inv_t = self.get_inverse_inertial_tensor_body_space();
        let orient = Mat3::from_quat(self.orientation);
        orient * inv_t * orient.transpose()
    }
    fn apply_impluse_linear(&mut self, impulse:Vec3) {
        if self.inv_mass == 0. {
            return;
        }
        self.linear_veclocity += impulse * self.inv_mass;
    }
    fn apply_impluse_angular(&mut self, impulse:Vec3) {
        if self.inv_mass == 0. {
            return;
        }
        self.angular_veclocity += self.get_inverse_inertial_tensor_world_space() * impulse;
        let max_angular_speed = 30.;
        if self.angular_veclocity.length_squared() > max_angular_speed * max_angular_speed {
            self.angular_veclocity = self.angular_veclocity.normalize() * max_angular_speed;
        }
    }
    fn apply_impulse(&mut self, impulse_point:Vec3, impulse:Vec3) {
        if self.inv_mass == 0. {
            return;
        }
        self.apply_impluse_linear(impulse);
        let r = impulse_point - self.get_centre_of_mass_world_space();
        self.apply_impluse_angular(r.cross(impulse));         
    }
    fn update(&mut self, dt_sec:f32) {
        self.position += self.linear_veclocity * dt_sec;

        let position_cm = self.get_centre_of_mass_world_space();
        let cm_to_pos = self.position - position_cm;
        let orient = Mat3::from_quat(self.orientation);
        let inertial_tensor = orient * self.shape.inertial_tensor() * orient.transpose();
        let alpha = inertial_tensor.inverse() *
            self.angular_veclocity.cross(inertial_tensor * self.angular_veclocity);
        self.angular_veclocity += alpha * dt_sec;

        let d_angle = self.angular_veclocity * dt_sec;
        let d_q = Quat::from_scaled_axis(d_angle);
        self.orientation = (d_q * self.orientation).normalize();

        self.position = position_cm + d_q.mul_vec3(cm_to_pos); 
    }
    fn radius(&self) -> f32 {
        self.shape.radius()
    }
    fn draw(&self, texture:Texture2D) {
        self.shape.draw(self.position, texture);
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
    let impulse_friction:Vec3;
    let pt_on_a = contact.pt_on_a_world_space;
    let pt_on_b = contact.pt_on_b_world_space;
    {
        let (a, b) = (&bodies[contact.body_a], &bodies[contact.body_b]);
        let elasticity = a.elasticity * b.elasticity;
        let inv_world_inertia_a = a.get_inverse_inertial_tensor_world_space();
        let inv_world_inertia_b = b.get_inverse_inertial_tensor_world_space();
        let n = contact.normal;
        let ra = pt_on_a - a.get_centre_of_mass_world_space();
        let rb = pt_on_b - b.get_centre_of_mass_world_space();
        let angular_ja = (inv_world_inertia_a * ra.cross(n)).cross(ra);
        let angular_jb = (inv_world_inertia_b * rb.cross(n)).cross(rb);
        let angular_factor = (angular_ja + angular_jb).dot(n);

        let vel_a = a.linear_veclocity + a.angular_veclocity.cross(ra);
        let vel_b = b.linear_veclocity + b.angular_veclocity.cross(rb);

        let vab = vel_a - vel_b;
        let impulse_j = (1. + elasticity) * vab.dot(contact.normal) /
            (a.inv_mass + b.inv_mass + angular_factor);
        vec_impulse_j = n * impulse_j;

        let friction = a.friction * b.friction;
        let vel_norm = n * n.dot(vab);
        let vel_tang = vab - vel_norm;
        let relative_vel_tang = vel_tang.normalize();

        let inertia_a = (inv_world_inertia_a * ra.cross(relative_vel_tang)).cross(ra);
        let inertia_b = (inv_world_inertia_b * rb.cross(relative_vel_tang)).cross(rb);
        let inv_inertia = (inertia_a + inertia_b).dot(relative_vel_tang);

        let reduced_mass = 1. / (a.inv_mass + b.inv_mass + inv_inertia);
        impulse_friction = vel_tang * reduced_mass * friction;
    }
    bodies[contact.body_a].apply_impulse(pt_on_a, vec_impulse_j * -1.);
    bodies[contact.body_b].apply_impulse(pt_on_b, vec_impulse_j * 1.);
    bodies[contact.body_a].apply_impulse(pt_on_a, impulse_friction * -1.);
    bodies[contact.body_b].apply_impulse(pt_on_b, impulse_friction * 1.);


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
                position: Vec3::new(0., 1., 0.),
                orientation: Quat::IDENTITY,
                linear_veclocity: Vec3::new(0., 0., 1.),
                angular_veclocity: Vec3::ZERO,
                inv_mass: 1.,
                elasticity: 0.,
                friction: 0.9,
                shape: Box::new(Sphere{radius:1., color:BLUE})
            }
        );
        bodies.push(
            Body {
                position: Vec3::new(0., -1000., 0.),
                orientation: Quat::IDENTITY,
                linear_veclocity: Vec3::ZERO,
                angular_veclocity: Vec3::ZERO,
                inv_mass: 0.,
                elasticity: 1.,
                friction: 0.,
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
            body.update(dt_sec);
        }
   }

    fn draw(&self, texture: Texture2D) {
        for body in &self.bodies {
            body.draw(texture);
        }
    }
}

#[macroquad::main("3D")]
async fn main() {
    let mut scene = Scene::new();

    let bytes: Vec<u8> = vec![255, 0, 0, 192, 0, 255, 0, 192, 0, 0, 255, 192, 255, 255, 255, 192];
    let texture = Texture2D::from_rgba8(2, 2, &bytes);

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

        scene.draw(texture);
        // Back to screen space, render some text

        set_default_camera();
        draw_text("WELCOME TO 3D WORLD", 10.0, 20.0, 30.0, BLACK);

        next_frame().await
    }
}