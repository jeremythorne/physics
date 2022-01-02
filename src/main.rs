use miniquad::*;

use glam::{vec3, Mat4, EulerRot};

mod main_pipe;
mod objects;
mod engine;

use main_pipe::MainPipe;
use objects::Object;
use physics::quad_verts;

struct PipeBind {
    pipe: Pipeline,
    bind: Bindings
}

struct Stage {
    main: MainPipe,
    main_bind: Bindings,
    verts_per_draw: u16,
    copy: PipeBind,
    scene: engine::Scene,
    rx: f32,
    ry: f32,
}

fn copy_pipe(ctx: &mut Context, tex:Texture) -> PipeBind {
    let (vertices, indices) = quad_verts();
    let vertex_buffer = Buffer::immutable(ctx, BufferType::VertexBuffer, &vertices);
    let index_buffer = Buffer::immutable(ctx, BufferType::IndexBuffer, &indices);

    let bind = Bindings {
        vertex_buffers: vec![vertex_buffer],
        index_buffer: index_buffer,
        images: vec![tex],
    };

    let shader = Shader::new(
        ctx,
        copy_to_screen_shader::VERTEX,
        copy_to_screen_shader::FRAGMENT,
        copy_to_screen_shader::meta(),
    )
    .unwrap();

    let pipe = Pipeline::new(
        ctx,
        &[BufferLayout::default()],
        &[
            VertexAttribute::new("pos", VertexFormat::Float2),
            VertexAttribute::new("uv", VertexFormat::Float2),
        ],
        shader,
    );

    PipeBind {
        pipe,
        bind
    }
}

impl Stage {
    pub fn new(ctx: &mut Context) -> Stage {
        let (bind, verts_per_draw) = objects::sphere_bindings(ctx);

        let main = MainPipe::new(ctx);
        let main_bind = bind.clone();

        let copy = copy_pipe(ctx, main.get_output());
        let scene = engine::Scene::new_cover();
 
        Stage {
            main,
            main_bind,
            verts_per_draw,
            copy,
            scene,
            rx: 0.,
            ry: 0.,
        }
    }
}

impl EventHandler for Stage {
    fn update(&mut self, _ctx: &mut Context) {
        self.scene.update(1.0 / 60.0);
    }

    fn resize_event(&mut self, ctx: &mut Context, width: f32, height: f32) {
        self.main.resize(ctx, width, height);
        self.copy.bind.images[0] = self.main.get_output();
    }

    fn draw(&mut self, ctx: &mut Context) {
        let (width, height) = ctx.screen_size();
        let proj = Mat4::perspective_rh_gl(45.0f32.to_radians(), width / height, 0.01, 100.0);
        let view = Mat4::look_at_rh(
            vec3(35.0, 20.0, 35.0),
            vec3(0.0, 0.0, 0.0),
            vec3(0.0, 1.0, 0.0),
        );
        let view_proj = proj * view;

        let proj = Mat4::perspective_rh_gl(60.0f32.to_radians(), 1.0, 10.0, 30.0);
        let light_view = Mat4::look_at_rh(
            vec3(10.0, 10.0, 10.0),
            vec3(0.0, 0.0, 0.0),
            vec3(0.0, 1.0, 0.0),
        );
        let light_view_proj = proj * light_view;
 
        //self.rx += 0.01;
        //self.ry += 0.01;
        let model = Mat4::from_euler(EulerRot::YXZ, self.ry, self.rx, 0.);

        //let (w, h) = ctx.screen_size();
       
        let matrices = self.scene.drawables();
        let mut objects = Vec::<Object>::new();
        for m in matrices {
            objects.push(Object {
                model: m,
                start: 0,
                end: self.verts_per_draw as i32,
            })
        }
        self.main.draw(ctx, &self.main_bind,
            &objects, 
            &model, &view_proj, &light_view_proj);

        let output = &self.copy;
        // and the post-processing-pass, rendering a quad, using the
        // previously rendered offscreen render-target as texture
        ctx.begin_default_pass(PassAction::Nothing);
        ctx.apply_pipeline(&output.pipe);
        ctx.apply_bindings(&output.bind);
        ctx.draw(0, 6, 1);
        ctx.end_render_pass();
        ctx.commit_frame();
    }
}

fn null_render() {
    let mut scene = engine::Scene::new_simple();
    let mut t = 0.;
    let dt = 1./60.;
    loop {
        scene.update(dt);
        t += dt;
        println!("time: {}s", t);
        for (i, body) in scene.bodies.iter().enumerate() {
            println!("{}: {} {} {} {} {} {} {}", i, body.position,
                body.orientation,
                body.linear_veclocity,
                body.angular_veclocity,
                body.inv_mass,
                body.elasticity,
                body.friction);
        }
    }
}

fn main() {
//    null_render();
    miniquad::start(conf::Conf::default(), |mut ctx| {
        UserData::owning(Stage::new(&mut ctx), ctx)
    });
}

mod copy_to_screen_shader {
    use miniquad::*;

    pub const VERTEX: &str = r#"#version 100
    attribute vec2 pos;
    attribute vec2 uv;

    varying lowp vec2 texcoord;

    void main() {
        gl_Position = vec4(pos, 0, 1);
        texcoord = uv;
    }
    "#;

    pub const FRAGMENT: &str = r#"#version 100
    precision lowp float;

    varying vec2 texcoord;

    uniform sampler2D tex;

    void main() {
        gl_FragColor = texture2D(tex, texcoord);
    }
    "#;

    pub fn meta() -> ShaderMeta {
        ShaderMeta {
            images: vec!["tex".to_string()],
            uniforms: UniformBlockLayout {
                uniforms: vec![],
            },
        }
    }
}

