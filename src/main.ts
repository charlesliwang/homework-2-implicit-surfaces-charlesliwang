import {vec2, vec3, vec4} from 'gl-matrix';
import * as Stats from 'stats-js';
import * as DAT from 'dat-gui';
import Square from './geometry/Square';
import Camera from './Camera';
import {setGL} from './globals';
import ShaderProgram, {Shader} from './rendering/gl/ShaderProgram';

// Define an object with application parameters and button callbacks
// This will be referred to by dat.GUI's functions that add GUI elements.
const controls = {
  'Duplicates' : false,
  // TODO: add any controls you want
};

let screenQuad: Square;
let time = 0;
let cycle = 300;

function main() {
  // Initial display for framerate
  const stats = Stats();
  stats.setMode(0);
  stats.domElement.style.position = 'absolute';
  stats.domElement.style.left = '0px';
  stats.domElement.style.top = '0px';
  document.body.appendChild(stats.domElement);

  // TODO: add any controls you need to the gui
  const gui = new DAT.GUI();
  gui.add(controls, 'Duplicates', false );
  // E.G. gui.add(controls, 'tesselations', 0, 8).step(1);

  // get canvas and webgl context
  const canvas = <HTMLCanvasElement> document.getElementById('canvas');

  function setSize(width: number, height: number) {
    canvas.width = width;
    canvas.height = height;
  }

  const gl = <WebGL2RenderingContext> canvas.getContext('webgl2');
  if (!gl) {
    alert('WebGL 2 not supported!');
  }
  // `setGL` is a function imported above which sets the value of `gl` in the `globals.ts` module.
  // Later, we can import `gl` from `globals.ts` to access it
  setGL(gl);

  screenQuad = new Square(vec3.fromValues(0, 0, 0));
  screenQuad.create();

  const camera = new Camera(vec3.fromValues(0, 0, 5), vec3.fromValues(0, 0, 0));

  gl.clearColor(0.0, 0.0, 0.0, 1);
  gl.disable(gl.DEPTH_TEST);

  const raymarchShader = new ShaderProgram([
    new Shader(gl.VERTEX_SHADER, require('./shaders/screenspace-vert.glsl')),
    new Shader(gl.FRAGMENT_SHADER, require('./shaders/raymarch-frag.glsl')),
  ]);

  // This function will be called every frame
  function tick() {
    camera.update();
    stats.begin();

    gl.viewport(0, 0, window.innerWidth, window.innerHeight);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

    // TODO: get / calculate relevant uniforms to send to shader here
    // TODO: send uniforms to shader

    // March!
    raymarchShader.setViewMatrix(camera.viewMatrix);
    time++;
    let dup = 0.0;
    if(controls.Duplicates) {
      dup = 1.0;
    }
    let v4 = vec4.fromValues(time/cycle,time%cycle,time,dup);
    raymarchShader.setTime(v4);
    raymarchShader.draw(screenQuad);

    // TODO: more shaders to layer / process the first one? (either via framebuffers or blending)

    stats.end();

    // Tell the browser to call `tick` again whenever it renders a new frame
    requestAnimationFrame(tick);
  }

  window.addEventListener('resize', function() {
    setSize(window.innerWidth, window.innerHeight);
    camera.setAspectRatio(window.innerWidth / window.innerHeight);
    camera.updateProjectionMatrix();
    let w = vec2.fromValues(window.innerWidth, window.innerHeight);
    raymarchShader.setWindow(w);
  }, false);

  setSize(window.innerWidth, window.innerHeight);
  camera.setAspectRatio(window.innerWidth / window.innerHeight);
  camera.updateProjectionMatrix();
  let w = vec2.fromValues(window.innerWidth, window.innerHeight);
  raymarchShader.setWindow(w);

  // Start the render loop
  tick();
}

main();
