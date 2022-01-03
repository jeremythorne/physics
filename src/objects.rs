use miniquad::*;

use glam::{Vec2, Vec3, Vec4, Mat4};

pub struct Drawable {
    pub model:Mat4,
    pub color:Vec4,
}

pub struct Object {
    pub drawable:Drawable,
    pub start:i32,
    pub end:i32
}

fn cube_verts() -> (&'static[f32], &'static[u16]) {
    #[rustfmt::skip]
    let vertices: &[f32] = &[
        /* pos               color                   uvs */
        -1.0, -1.0, -1.0,    1.0, 0.5, 0.5, 1.0,     0.0, 0.0,
        1.0, -1.0, -1.0,    1.0, 0.5, 0.5, 1.0,     1.0, 0.0,
        1.0,  1.0, -1.0,    1.0, 0.5, 0.5, 1.0,     1.0, 1.0,
        -1.0,  1.0, -1.0,    1.0, 0.5, 0.5, 1.0,     0.0, 1.0,

        -1.0, -1.0,  1.0,    0.5, 1.0, 0.5, 1.0,     0.0, 0.0,
        1.0, -1.0,  1.0,    0.5, 1.0, 0.5, 1.0,     1.0, 0.0,
        1.0,  1.0,  1.0,    0.5, 1.0, 0.5, 1.0,     1.0, 1.0,
        -1.0,  1.0,  1.0,    0.5, 1.0, 0.5, 1.0,     0.0, 1.0,

        -1.0, -1.0, -1.0,    0.5, 0.5, 1.0, 1.0,     0.0, 0.0,
        -1.0,  1.0, -1.0,    0.5, 0.5, 1.0, 1.0,     1.0, 0.0,
        -1.0,  1.0,  1.0,    0.5, 0.5, 1.0, 1.0,     1.0, 1.0,
        -1.0, -1.0,  1.0,    0.5, 0.5, 1.0, 1.0,     0.0, 1.0,

        1.0, -1.0, -1.0,    1.0, 0.5, 0.0, 1.0,     0.0, 0.0,
        1.0,  1.0, -1.0,    1.0, 0.5, 0.0, 1.0,     1.0, 0.0,
        1.0,  1.0,  1.0,    1.0, 0.5, 0.0, 1.0,     1.0, 1.0,
        1.0, -1.0,  1.0,    1.0, 0.5, 0.0, 1.0,     0.0, 1.0,

        -1.0, -1.0, -1.0,    0.0, 0.5, 1.0, 1.0,     0.0, 0.0,
        -1.0, -1.0,  1.0,    0.0, 0.5, 1.0, 1.0,     1.0, 0.0,
        1.0, -1.0,  1.0,    0.0, 0.5, 1.0, 1.0,     1.0, 1.0,
        1.0, -1.0, -1.0,    0.0, 0.5, 1.0, 1.0,     0.0, 1.0,

        -1.0,  1.0, -1.0,    1.0, 0.0, 0.5, 1.0,     0.0, 0.0,
        -1.0,  1.0,  1.0,    1.0, 0.0, 0.5, 1.0,     1.0, 0.0,
        1.0,  1.0,  1.0,    1.0, 0.0, 0.5, 1.0,     1.0, 1.0,
        1.0,  1.0, -1.0,    1.0, 0.0, 0.5, 1.0,     0.0, 1.0
    ];

    #[rustfmt::skip]
    let indices: &[u16] = &[
        0, 1, 2,  0, 2, 3,
        6, 5, 4,  7, 6, 4,
        8, 9, 10,  8, 10, 11,
        14, 13, 12,  15, 14, 12,
        16, 17, 18,  16, 18, 19,
        22, 21, 20,  23, 22, 20
    ];

    (vertices, indices)
}

fn bindings(ctx: &mut Context, verts:(&[f32],&[u16])) -> Bindings {
    let (vertices, indices) = verts;
    let vertex_buffer = Buffer::immutable(ctx, BufferType::VertexBuffer, &vertices);
    let index_buffer = Buffer::immutable(ctx, BufferType::IndexBuffer, &indices);

    Bindings {
        vertex_buffers: vec![vertex_buffer],
        index_buffer: index_buffer,
        images: vec![],
    }
}

pub fn cube_bindings(ctx: &mut Context) -> Bindings {
    bindings(ctx, cube_verts())
}

fn octahedral_verts() -> (&'static[f32], &'static[u16]) {
    #[rustfmt::skip]
    let vertices: &[f32] = &[
        /* pos               color                   uvs */
        0.0,  1.0,  0.0,    1.0, 0.5, 0.5, 1.0,     0.125, 0.0,
        -1.0, 0.0,  0.0,    1.0, 0.5, 0.5, 1.0,     0.0,   0.5,
        0.0,  0.0,  1.0,    1.0, 0.5, 0.5, 1.0,     0.25,  0.5,

        0.0,  1.0,  0.0,    0.5, 1.0, 0.5, 1.0,     0.125, 0.0,
        0.0,  0.0,  1.0,    0.5, 1.0, 0.5, 1.0,     0.25,  0.5,
        1.0,  0.0,  0.0,    0.5, 1.0, 0.5, 1.0,     0.5,   0.5,

        0.0,  1.0,  0.0,    0.5, 0.5, 1.0, 1.0,     0.125, 0.0,
        1.0,  0.0,  0.0,    0.5, 0.5, 1.0, 1.0,     0.5,   0.5,
        0.0,  0.0, -1.0,    0.5, 0.5, 1.0, 1.0,     0.75,  0.5,

        0.0,  1.0,  0.0,    1.0, 0.0, 0.0, 1.0,     0.125, 0.0,
        0.0,  0.0, -1.0,    1.0, 0.0, 0.0, 1.0,     0.75,  0.5,
       -1.0,  0.0,  0.0,    1.0, 0.0, 0.0, 1.0,     1.0,   0.5,

        0.0, -1.0,  0.0,    0.0, 1.0, 0.0, 1.0,     0.125, 1.0,
       -1.0,  0.0,  0.0,    0.0, 1.0, 0.0, 1.0,     0.0,   0.5,
        0.0,  0.0,  1.0,    0.0, 1.0, 0.0, 1.0,     0.25,  0.5,

        0.0, -1.0,  0.0,    0.0, 0.0, 1.0, 1.0,     0.125, 1.0,
        0.0,  0.0,  1.0,    0.0, 0.0, 1.0, 1.0,     0.25,  0.5,
        1.0,  0.0,  0.0,    0.0, 0.0, 1.0, 1.0,     0.5,   0.5,

        0.0, -1.0,  0.0,    1.0, 1.0, 0.0, 1.0,     0.125, 1.0,
        1.0,  0.0,  0.0,    1.0, 1.0, 0.0, 1.0,     0.5,   0.5,
        0.0,  0.0, -1.0,    1.0, 1.0, 0.0, 1.0,     0.75,  0.5,

        0.0, -1.0,  0.0,    0.0, 1.0, 1.0, 1.0,     0.125, 1.0,
        0.0,  0.0, -1.0,    0.0, 1.0, 1.0, 1.0,     0.75,  0.5,
       -1.0,  0.0,  0.0,    0.0, 1.0, 1.0, 1.0,     1.0,   0.5,
    ];

    #[rustfmt::skip]
    let indices: &[u16] = &[
        0,  1,  2,
        3,  4,  5,
        6,  7,  8,
        9, 10, 11,
       14, 13, 12,
       17, 16, 15,
       20, 19, 18,
       23, 22, 21 
    ];

    (vertices, indices)
}

pub fn octahedral_bindings(ctx: &mut Context) -> Bindings {
    let (vertices, indices) = octahedral_verts();
    bindings(ctx, (&vertices, &indices))
}


//from https://sites.google.com/site/dlampetest/python/triangulating-a-sphere-recursively
// Subdivide each triangle in the old approximation and normalize
//  the new points thus generated to lie on the surface of the unit
//  sphere.
// Each input triangle with vertices labelled [0,1,2] as shown
//  below will be turned into four new triangles:
//
//            Make new points     
//                 a = (0+2)/2
//                 b = (0+1)/2     
//                 c = (1+2)/2
//        1     
//       /\        Normalize a, b, c     
//      /  \     
//    b/____\ c    Construct new new triangles     
//    /\    /\       t1 [0,b,a]     
//   /  \  /  \      t2 [b,1,c]     
//  /____\/____\     t3 [a,b,c]     
// 0      a     2    t4 [a,c,2]

fn average_vert(a:&[f32], b:&[f32]) -> Vec<f32> {
    let pos_a = Vec3::from_slice(&a[0..3]);
    let col_a = Vec4::from_slice(&a[3..7]);
    let uv_a = Vec2::from_slice(&a[7..9]);
    let pos_b = Vec3::from_slice(&b[0..3]);
    let col_b = Vec4::from_slice(&b[3..7]);
    let uv_b = Vec2::from_slice(&b[7..9]);

    let pos = ((pos_a + pos_b) / 2.).normalize();
    let col = (col_a + col_b) / 2.;
    let uv = (uv_a + uv_b) / 2.;
    vec![pos.x, pos.y, pos.z, col.x, col.y, col.z, col.w, uv.x, uv.y]
}

fn divide_one(vertices:&[f32], indices:&[u16]) -> (Vec<f32>, Vec<u16>, u16) {
    let stride = 9;
    let i0 = stride * indices[0] as usize; 
    let i1 = stride * indices[1] as usize;  
    let i2 = stride * indices[2] as usize;
    let v0 = &vertices[i0..i0 + stride];
    let v1 = &vertices[i1..i1 + stride];
    let v2 = &vertices[i2..i2 + stride];
    let a = average_vert(v0, v2);
    let b = average_vert(v0, v1);
    let c = average_vert(v1, v2);
    let mut vertices_out = Vec::<f32>::new();
    let max_verts = 12;
    let indices_out: Vec<u16> = (0..max_verts).collect();

    let mut push = |v: &[&[f32]]| {
        for vv in v {
            vertices_out.extend_from_slice(vv)
        }
    };
    push(&[v0, &b, &a]);
    push(&[&b, v1, &c]);
    push(&[&a, &b, &c]);
    push(&[&a, &c, v2]);

    (vertices_out, indices_out, max_verts)
}

fn divide_all(vertices:&[f32], indices:&[u16]) -> (Vec<f32>, Vec<u16>, u16) {
    let mut vertices_out = Vec::<f32>::new();
    let mut indices_out = Vec::<u16>::new();
    let mut max_verts = 0;
    let n = indices.len() / 3;
    for i in 0..n {
        let (vo, io, mv) = divide_one(vertices, &indices[i * 3..i * 3 + 3]);
        let no = indices_out.len();
        vertices_out.extend_from_slice(&vo);
        let mut iio = Vec::<u16>::new();
        for i in io {
            iio.push(i + (no as u16));
        }
        indices_out.extend_from_slice(&iio);
        max_verts += mv;
    }
    (vertices_out, indices_out, max_verts)
} 

fn sphere_verts() -> (Vec<f32>,Vec<u16>, u16) {
    let (vertices, indices) = octahedral_verts();
    let mut vertices_out = Vec::<f32>::from(vertices);
    let mut indices_out = Vec::<u16>::from(indices);
    let mut max_verts = 0;
    for _ in 0..3 {
        let (vo, io, mv) = divide_all(&vertices_out, &indices_out);
        vertices_out = vo;
        indices_out = io;
        max_verts = mv;
    }
    (vertices_out, indices_out, max_verts)
}

pub fn sphere_bindings(ctx: &mut Context) -> (Bindings, u16) {
    let (vertices, indices, max_verts) = sphere_verts();
    (bindings(ctx, (&vertices, &indices)), max_verts)
}

