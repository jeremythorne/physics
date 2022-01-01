pub fn quad_verts() -> (&'static[f32], &'static[u16]) {
    #[rustfmt::skip]
    let vertices: &[f32] = &[
        /* pos         uvs */
        -1.0, -1.0,    0.0, 0.0,
        1.0, -1.0,    1.0, 0.0,
        1.0,  1.0,    1.0, 1.0,
        -1.0,  1.0,    0.0, 1.0,
    ];
    let indices: &[u16] = &[0, 1, 2, 0, 2, 3];
    (vertices, indices)
}