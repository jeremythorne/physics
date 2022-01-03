use glam::{Vec3, Vec4, Mat3, Mat4, Quat};
use std::cmp::Ordering;
use std::any::Any;
use std::fmt::Debug;
use crate::objects::Drawable;

struct Sphere {
    radius: f32,
    color: Vec4,
}

#[derive(PartialEq)]
enum Shape {
    Sphere
}

trait Shaped {
    fn shape_type(&self) -> Shape;
    fn color(&self) -> Vec4;
    fn get_bounds(&self) -> Bounds;
    fn get_bounds_pos(&self, pos: Vec3, orient: Quat) -> Bounds;
    fn get_centre_of_mass(&self) -> Vec3;
    fn inertial_tensor(&self) -> Mat3;
    fn radius(&self) -> f32;
    fn as_any(&self) -> &dyn Any;
}

impl Shaped for Sphere {
    fn shape_type(&self) -> Shape {
        Shape::Sphere
    }

    fn color(&self) -> Vec4 {
        self.color
    }

    fn get_bounds(&self) -> Bounds {
        Bounds {
            mins: Vec3::splat(-self.radius()),
            maxs: Vec3::splat(self.radius())
        }
    }

    fn get_bounds_pos(&self, pos: Vec3, _: Quat) -> Bounds {
        Bounds {
            mins: Vec3::splat(-self.radius()) + pos,
            maxs: Vec3::splat(self.radius()) + pos
        }
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
}

pub struct Body {
    pub position: Vec3,
    pub orientation: Quat,
    pub linear_veclocity: Vec3,
    pub angular_veclocity: Vec3,
    pub inv_mass: f32,
    pub elasticity: f32,
    pub friction: f32,
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
}

pub struct Scene {
    pub bodies: Vec<Body>
}

struct Bounds {
    mins: Vec3,
    maxs: Vec3
}

impl Bounds {
    fn new() -> Bounds {
        Bounds {
            mins: Vec3::splat(1e6),
            maxs: Vec3::splat(-1e6)
        }
    }
    fn does_intersect(self:&Bounds, other:&Bounds) -> bool {
        if self.maxs.x < other.mins.x || self.maxs.y < other.mins.y || self.maxs.z < other.mins.z {
            return false;
        }
        if other.maxs.x < self.mins.x || other.maxs.y < self.mins.y || other.maxs.z < self.mins.z {
            return false;
        }
        true
    }
    fn expand_vec(&mut self, bounds: &[Bounds]) {
        for b in bounds {
            self.expand_bound(b);
        }
    }
    fn expand_bound(&mut self, other: &Bounds) {
        self.expand(other.mins);
        self.expand(other.maxs);
    }
    fn expand(&mut self, v: Vec3) {
        self.mins = self.mins.min(v);
        self.maxs = self.maxs.max(v);
    }
}

struct PseudoBody {
    id: usize,
    value: f32,
    is_min: bool
}

fn sort_bodies_bounds(bodies:&[Body], dt_sec:f32) -> Vec<PseudoBody> {
    let axis = Vec3::new(1., 1., 1.).normalize();
    let mut sorted = Vec::<PseudoBody>::new();

    for i in 0..bodies.len() {
        let body = &bodies[i];
        let mut bounds = body.shape.get_bounds_pos(body.position, body.orientation);
        bounds.expand(bounds.mins + body.linear_veclocity * dt_sec);
        bounds.expand(bounds.maxs + body.linear_veclocity * dt_sec);
        let epsilon = 0.01;
        bounds.expand(bounds.mins + Vec3::splat(-1.) * epsilon);
        bounds.expand(bounds.maxs + Vec3::splat(1.) * epsilon);

        sorted.push(PseudoBody{
            id: i,
            value: axis.dot(bounds.mins),
            is_min: true
        });
        sorted.push(PseudoBody{
            id: i,
            value: axis.dot(bounds.maxs),
            is_min: false
        });
    }
    sorted.sort_by(|a, b|
        a.value.partial_cmp(&b.value).unwrap_or(Ordering::Equal));
    sorted
}

#[derive(PartialEq)]
struct CollisionPair {
    a: usize,
    b: usize
}

fn build_pairs(sorted_bodies:&[PseudoBody]) -> Vec<CollisionPair> {
    let mut collision_pairs = Vec::<CollisionPair>::new();
    for i in 0..sorted_bodies.len() {
        let a = &sorted_bodies[i];
        if !a.is_min {
            continue;
        }
        for j in i + 1..sorted_bodies.len() {
            let b = &sorted_bodies[j];
            if b.id == a.id {
                // done with a
                break;
            }
            if !b.is_min {
                continue;
            }
            collision_pairs.push(CollisionPair{
                a: a.id,
                b: b.id
            })
        }
    }
    collision_pairs
}

fn sweep_and_prune_1d(bodies:&[Body], dt_sec:f32) -> Vec<CollisionPair> {
    let sorted_bodies = sort_bodies_bounds(bodies, dt_sec);
    build_pairs(&sorted_bodies)
}

fn broad_phase(bodies:&[Body], dt_sec:f32) -> Vec<CollisionPair> {
    sweep_and_prune_1d(bodies, dt_sec)
}

#[derive(Debug)]
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

fn ray_sphere(ray_start:Vec3, ray_dir:Vec3, sphere_centre:Vec3, radius:f32) -> 
        Option<(f32, f32)> {
    let m = sphere_centre - ray_start;
    let a = ray_dir.dot(ray_dir);
    let b = m.dot(ray_dir);
    let c = m.dot(m) - radius * radius;
    let delta = b * b - a * c;
    let inv_a = 1. / a;
    if delta < 0. {
        return None;
    }
    let delta_root = delta.sqrt();
    let t1 = inv_a * (b - delta_root);
    let t2 = inv_a * (b + delta_root);
    Some((t1, t2))
}

fn sphere_sphere_dynamic(shape_a:&dyn Shaped, shape_b:&dyn Shaped, pos_a:Vec3, pos_b:Vec3, 
    vel_a:Vec3, vel_b:Vec3, dt:f32) ->
        Option<(Vec3, Vec3, f32)> {

    let relative_velocity = vel_a - vel_b;
    let ray_dir = relative_velocity * dt;
    let mut t0 = 0.;
    let mut t1 = 0.;
    if ray_dir.length_squared() < 0.001 * 0.001 {
        let ab = pos_b - pos_a;
        let radius = shape_a.radius() + shape_b.radius() + 0.001;
        if ab.length_squared() > radius * radius {
            return None;
        }
    } else {
        if let Some(x) = ray_sphere(pos_a, ray_dir, pos_b, shape_a.radius() + shape_b.radius()) {
            t0 = x.0;
            t1 = x.1;
        } else {
            return None;
        }
    }
    t0 *= dt;
    t1 *= dt;
    if t1 < 0. {
        return None;
    }
    let time_of_impact = t0.max(0.);
    if time_of_impact > dt {
        return None;
    }

    let new_pos_a = pos_a + vel_a * time_of_impact;
    let new_pos_b = pos_b + vel_b * time_of_impact;
    let ab = (new_pos_b - new_pos_a).normalize();

    let pt_on_a = new_pos_a + ab * shape_a.radius();
    let pt_on_b = new_pos_b - ab * shape_b.radius();

    Some((pt_on_a, pt_on_b, time_of_impact))
}

fn intersect(bodies: &mut [Body], i: usize, j: usize, dt:f32) -> Option<Contact> {
    let (a, b) = (&bodies[i], &bodies[j]);
    if a.shape.shape_type() == Shape::Sphere &&
            b.shape.shape_type() == Shape::Sphere {
        if let Some(x) = sphere_sphere_dynamic(
                &*a.shape, &*b.shape, a.position, b.position,
                a.linear_veclocity,
                b.linear_veclocity, dt) {
            let pt_on_a_world_space = x.0;
            let pt_on_b_world_space = x.1;
            let time_of_impact = x.2;
            let ab = b.position - a.position;
            let separation_distance = ab.length() - (a.radius() + b.radius());
            // step forward to calculate local coords
            bodies[i].update(time_of_impact);
            bodies[j].update(time_of_impact);
            let (a, b) = (&bodies[i], &bodies[j]);
            let pt_on_a_local_space = a.world_space_to_body_space(pt_on_a_world_space);
            let pt_on_b_local_space = b.world_space_to_body_space(pt_on_b_world_space);
            let normal = (a.position - b.position).normalize();
            // step backward
            bodies[i].update(-time_of_impact);
            bodies[j].update(-time_of_impact);
            return Some(Contact {
                pt_on_a_world_space,
                pt_on_b_world_space,
                pt_on_a_local_space,
                pt_on_b_local_space,
                normal,
                separation_distance,
                time_of_impact,
                body_a: i,
                body_b: j
            });
        } else {
            return None;
        }
    }
    None
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
        let relative_vel_tang = vel_tang.normalize_or_zero();

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


    if contact.time_of_impact == 0. {
        // move colliders to just outside each other but keeping combined centre
        // of mass constant
        let (a, b) = (&bodies[contact.body_a], &bodies[contact.body_b]);
        let ta = a.inv_mass / (a.inv_mass + b.inv_mass);
        let tb = b.inv_mass / (a.inv_mass + b.inv_mass);
        let ds = contact.pt_on_b_world_space - contact.pt_on_a_world_space;
        bodies[contact.body_a].position += ds * ta;
        bodies[contact.body_b].position -= ds * tb;
    }
}

impl Scene {
    pub fn new_simple() -> Scene {
        let mut bodies = Vec::<Body>::new();
        // dynamic bodies
        let radius = 2.0;
        let x = 0.;
        let z = 0.;
        let y = 10.;
        let color = Vec4::new(1.0, 0.2, 0.2, 1.0);
        bodies.push(
            Body {
                position: Vec3::new(z, y, x),
                orientation: Quat::IDENTITY,
                linear_veclocity: Vec3::ZERO,
                angular_veclocity: Vec3::ZERO,
                inv_mass: 1.0,
                elasticity: 1.,
                friction: 0.5,
                shape: Box::new(Sphere{radius, color})
            }
        );


        // static "floor"
        let radius = 100.;
        let x = 0.;
        let z = 0.;
        let y = -radius;
        let color = Vec4::new(0.2, 1.0, 0.2, 1.0);
        bodies.push(
            Body {
                position: Vec3::new(z, y, x),
                orientation: Quat::IDENTITY,
                linear_veclocity: Vec3::ZERO,
                angular_veclocity: Vec3::ZERO,
                inv_mass: 0.,
                elasticity: 1.,
                friction: 0.5,
                shape: Box::new(Sphere{radius, color})
            }
        );
        Scene {
            bodies
        }
    }

    pub fn new_cover() -> Scene {
        let mut bodies = Vec::<Body>::new();
        // dynamic bodies
        for i in 0..6 {
            for j in 0..6 {
                for k in 0..6 {
                    let radius = 0.5;
                    let x = ((i - 1) as f32) * radius * 1.5;
                    let z = ((j - 1) as f32) * radius * 1.5;
                    let y = 10. + radius * 1.5 * k as f32;
                    let color = Vec4::new(1.0, 0.2, 0.2, 1.0);
                    bodies.push(
                    Body {
                            position: Vec3::new(z, y, x),
                            orientation: Quat::IDENTITY,
                            linear_veclocity: Vec3::ZERO,
                            angular_veclocity: Vec3::ZERO,
                            inv_mass: 1.,
                            elasticity: 0.5,
                            friction: 0.5,
                            shape: Box::new(Sphere{radius, color})
                        }
                    );
                }
            }
        }

        // static "floor"
        for i in 0..3 {
            for j in 0..3 {
                let radius = 80.;
                let x = ((i - 1) as f32) * radius * 0.25;
                let z = ((j - 1) as f32) * radius * 0.25;
                let y = -radius;
                let color = Vec4::new(0.2, 1.0, 0.2, 1.0);
                bodies.push(
                Body {
                        position: Vec3::new(z, y, x),
                        orientation: Quat::IDENTITY,
                        linear_veclocity: Vec3::ZERO,
                        angular_veclocity: Vec3::ZERO,
                        inv_mass: 0.,
                        elasticity: 0.99,
                        friction: 0.5,
                        shape: Box::new(Sphere{radius, color})
                    }
                );
            }
        }

        Scene {
            bodies
        }
    }

    pub fn new_many() -> Scene {
        let mut bodies = Vec::<Body>::new();
        // dynamic bodies
        for i in 0..6 {
            for j in 0..6 {
                let radius = 0.5;
                let x = ((i - 1) as f32) * radius * 1.5;
                let z = ((j - 1) as f32) * radius * 1.5;
                let y = 10.;
                let color = Vec4::new(1.0, 0.2, 0.2, 1.0);
                bodies.push(
                Body {
                        position: Vec3::new(z, y, x),
                        orientation: Quat::IDENTITY,
                        linear_veclocity: Vec3::ZERO,
                        angular_veclocity: Vec3::ZERO,
                        inv_mass: 1.,
                        elasticity: 0.5,
                        friction: 0.5,
                        shape: Box::new(Sphere{radius, color})
                    }
                );
            }
        }

        // static "floor"
        for i in 0..3 {
            for j in 0..3 {
                let radius = 80.;
                let x = ((i - 1) as f32) * radius * 0.25;
                let z = ((j - 1) as f32) * radius * 0.25;
                let y = -radius;
                let color = Vec4::new(0.2, 1.0, 0.2, 1.0);
                bodies.push(
                Body {
                        position: Vec3::new(z, y, x),
                        orientation: Quat::IDENTITY,
                        linear_veclocity: Vec3::ZERO,
                        angular_veclocity: Vec3::ZERO,
                        inv_mass: 0.,
                        elasticity: 0.99,
                        friction: 0.5,
                        shape: Box::new(Sphere{radius, color})
                    }
                );
            }
        }

        Scene {
            bodies
        }
    }

    pub fn update(&mut self, dt_sec: f32) {
        // gravity
        for body in &mut self.bodies {
            let mass = 1. / body.inv_mass;
            let impulse_gravity = Vec3::new(0., -10., 0.) * mass * dt_sec;
            body.apply_impluse_linear(impulse_gravity);
        }
        //broadphase
        let collision_pairs = broad_phase(&self.bodies, dt_sec);

        //narrow phase (perform actual collision detection)
        let mut contacts = Vec::<Contact>::new();
        for pair in collision_pairs {
            // skip body pairs with infinite mass
            let (a, b) = (&self.bodies[pair.a], &self.bodies[pair.b]);
            if a.inv_mass == 0. && b.inv_mass == 0. {
                continue
            }

            if let Some(contact) = intersect(&mut self.bodies[..], pair.a, pair.b, dt_sec) {
                contacts.push(contact);
            }
        }

        // sort by increasing time of impact
        contacts.sort_by(|a, b| a.time_of_impact
            .partial_cmp(&b.time_of_impact).unwrap());

        // apply ballistic impulses
        let mut accumulated_time = 0.;
        for contact in &contacts {
            let dt = contact.time_of_impact - accumulated_time;
            let (a, b) = (&self.bodies[contact.body_a], &self.bodies[contact.body_b]);
            if a.inv_mass == 0. && b.inv_mass == 0. {
                continue
            }
            // position update
            for body in &mut self.bodies {
                body.update(dt);
            }
            resolve_contact(&mut self.bodies[..], &contact);
            accumulated_time += dt;
        }

        // position update for time remaining in this frame
        let time_remaining = dt_sec - accumulated_time;
        for body in &mut self.bodies {
            body.update(time_remaining);
        }
   }

    pub fn drawables(&self) -> Vec<Drawable> {
        let mut objects = Vec::<Drawable>::new();
        for body in &self.bodies {
            let rot = Mat4::from_quat(body.orientation);
            let scale = Mat4::from_scale(Vec3::splat(body.shape.radius()));
            let trans = Mat4::from_translation(body.position);
            objects.push(Drawable {
                model: trans * rot * scale,
                color: body.shape.color()
            });
        }
        objects
    }
}
