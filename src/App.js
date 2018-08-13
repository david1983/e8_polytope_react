import React, { Component } from 'react';
import './App.css';

import { observable } from "mobx"
import { observer } from "mobx-react"

// Coordinates of the roots
var roots = new Array();

var canvas0, canvas1, ctx0, ctx1, candiv
const width = window.innerWidth;
const height = window.innerHeight;
var zoom = 100
var translateX = width / 2
var translateY = height / 2


function createRoots() {
  function rootFirstKind(i0, i1, s0, s1) {
    var rt = [0, 0, 0, 0, 0, 0, 0, 0];
    rt[i0] = (s0 ? -2 : 2);
    rt[i1] = (s1 ? -2 : 2);
    return rt;
  }
  function rootSecondKind(s0, s1, s2, s3, s4, s5, s6) {
    var s7 = s0 ^ s1 ^ s2 ^ s3 ^ s4 ^ s5 ^ s6;
    return [(s0 ? -1 : 1), (s1 ? -1 : 1), (s2 ? -1 : 1), (s3 ? -1 : 1),
    (s4 ? -1 : 1), (s5 ? -1 : 1), (s6 ? -1 : 1), (s7 ? -1 : 1)];
  }
  for (var i0 = 0; i0 < 8; i0++) {
    for (var i1 = i0 + 1; i1 < 8; i1++) {
      roots.push(rootFirstKind(i0, i1, false, false));
      roots.push(rootFirstKind(i0, i1, false, true));
      roots.push(rootFirstKind(i0, i1, true, false));
      roots.push(rootFirstKind(i0, i1, true, true));
    }
  }
  for (var i = 0; i < 128; i++) {
    roots.push(rootSecondKind(!!(i & 1), !!(i & 2), !!(i & 4), !!(i & 8),
      !!(i & 16), !!(i & 32), !!(i & 64)));
  }
  roots.sort(function (a, b) { // Lexicographic ordering
    for (var k = 0; k < 8; k++) {
      if (a[k] < b[k])
        return -1;
      else if (a[k] > b[k])
        return 1;
    }
    return 0;
  });
}

function sqnorm(a) {
  var ret = 0;
  for (var i = 0; i < 8; i++) {
    var d = a[i];
    ret += d * d;
  }
  return ret;
}

function dotprod(a, b) {
  var ret = 0;
  for (var i = 0; i < 8; i++) {
    ret += a[i] * b[i];
  }
  return ret;
}

function sqdist(a, b) {
  var ret = 0;
  for (var i = 0; i < 8; i++) {
    var d = a[i] - b[i];
    ret += d * d;
  }
  return ret;
}

// List of roots connected to each given root by an edge:
var adjacent = new Array();

function createAdjacent() {
  for (var a = 0; a < roots.length; a++) {
    adjacent[a] = new Array();
    for (var b = 0; b < roots.length; b++) {
      if (sqdist(roots[a], roots[b]) == 8)
        adjacent[a].push(b);
    }
  }
}

var gaussianStore = null;
function gaussian() {  // Generate a Gaussian variable by Box-Muller.
  if (gaussianStore == null) {
    var u0 = Math.random();
    var u1 = Math.random();
    gaussianStore = Math.sqrt(-2 * Math.log(u0)) * Math.cos(2 * Math.PI * u1);
    return Math.sqrt(-2 * Math.log(u0)) * Math.sin(2 * Math.PI * u1);
  } else {
    var ret = gaussianStore;
    gaussianStore = null;
    return ret;
  }
}

// The two 8-vectors determining the projection to a plane:
// these should be normed and orthogonal.
var projMatrix = observable(new Array());

function gramSchmidt() { // Make projMatrix normed and orthogonal.
  var ortho = [[], [0]];
  for (var k = 0; k < ortho.length; k++) {
    var d;
    for (var l = 0; l < ortho[k].length; l++) {
      var k2 = ortho[k][l];
      d = dotprod(projMatrix[k], projMatrix[k2]);
      for (var i = 0; i < 8; i++)
        projMatrix[k][i] -= d * projMatrix[k2][i];
    }
    d = Math.sqrt(sqnorm(projMatrix[k]));
    for (var i = 0; i < 8; i++)
      projMatrix[k][i] /= d;
  }
}

// The rotate8a, rotate8b, rotate30, etc., variables are lists of flips
// (wrt fundamental roots) to compose to create various Weyl group elements.

// This class has order 8 and centralizer of order 192, and acts on the root by 30 cycles of length 8.
var rotate8a = [150, 140, 132, 126, 122, 120, 162, 120, 122, 121, 126, 122, 120, 132, 126, 122, 121, 140, 132, 126, 122, 120, 162, 120, 122, 121, 126, 122, 120, 132, 126, 122, 121, 140, 132, 126, 122, 120, 150, 140, 132, 126, 122, 121, 162, 120];
// This class has order 8 and centralizer of order 64, and acts on the root by 27 cycles of length 8, 5 of length 4, 1 of length 2 and 2 fixed points.
var rotate8b = [162, 120, 122, 121, 126, 122, 120, 132, 126, 122, 121, 140, 132, 126, 122, 120, 162, 120, 122, 121, 126, 122, 120, 132, 126, 122, 121, 140, 132, 126, 122, 120, 150, 140, 132, 126, 122, 162, 120];

function chooseProjection8() {
  projMatrix[0] = new Array(8);
  projMatrix[1] = new Array(8);
  for (var k = 0; k < 8; k++) {
    projMatrix[0][k] = Math.cos((2 * k + 1) * Math.PI / 16.);
    projMatrix[1][k] = Math.sin((2 * k + 1) * Math.PI / 16.);
  }
  gramSchmidt();  // (Actually just divides by 2...)
}

function chooseProjectionAlt() {
  projMatrix[0] = [73. / 105., 67. / 105., 53. / 105., 17. / 105., -17. / 105., -53. / 105., -67. / 105., -73. / 105.];
  projMatrix[1] = [17. / 105., 53. / 105., 67. / 105., 73. / 105., 73. / 105., 67. / 105., 53. / 105., 17. / 105.];
  gramSchmidt();
}

function chooseProjectionSquares() {
  projMatrix[0] = [14, 8, 7, -1, -2, -4, -11, -13];
  projMatrix[1] = [2, 4, 11, 13, 14, 8, 7, -1];
  gramSchmidt();
}

// This class has order 30 and centralizer of order 30, and acts on the root by 8 cycles of length 30.
var rotate30 = [120, 121, 122, 120, 126, 122, 121, 132, 140, 150, 140, 162, 120, 122, 121, 126, 122, 120, 132, 126, 122, 121, 140, 132, 126, 122, 120, 150, 140, 132, 126, 122, 162, 120, 122, 121, 126, 122, 120, 132, 126, 122, 121, 140, 132, 126, 122, 120, 150, 140, 132, 126, 162, 120, 122, 121, 126, 122];

function chooseProjection30() {
  projMatrix[0] = [0.438217070641, 0.205187681291, 0.36459828198, 0.0124511903657, -0.0124511903657, -0.36459828198, -0.205187681291, -0.67645247517];
  projMatrix[1] = [-0.118465163028, 0.404927414852, 0.581970822973, 0.264896157496, 0.501826483552, 0.345040496917, 0.167997088796, 0.118465163028];
  gramSchmidt();  // (Actually unnecessary)
}

// This class has order 24 and centralizer of order 24, and acts on the root by 10 cycles of length 24.
var rotate24 = [121, 122, 120, 126, 122, 132, 150, 140, 132, 126, 122, 162, 120, 122, 121, 126, 122, 120, 132, 126, 122, 121, 140, 132, 126, 122, 120, 150, 140, 132, 126, 122, 162, 120, 122, 121, 126, 122, 120, 132, 126, 122, 140, 132, 126, 150, 162, 120, 122, 121];

function chooseProjection24() {
  projMatrix[0] = [0.622590719424, 0.197882519146, 0.382279671629, -0.0424285286489, -0.0819656231869, -0.257885519793, -0.498196567588, -0.322276670982];
  projMatrix[1] = [0.0819656231869, 0.257885519793, 0.498196567588, 0.322276670982, 0.622590719424, 0.197882519146, 0.382279671629, -0.0424285286489];
  gramSchmidt();  // (Actually unnecessary)
}

// This class has order 20 and centralizer of order 20, and acts on the root by 12 cycles of length 20.
var rotate20 = [120, 121, 122, 120, 121, 126, 122, 132, 126, 122, 140, 162, 120, 122, 121, 126, 122, 120, 132, 126, 122, 121, 140, 132, 126, 122, 120, 150, 140, 132, 126, 122, 162, 120, 122, 121, 126, 122, 120, 132, 126, 122, 121, 140, 132, 126, 162, 120, 122, 121];

function chooseProjection20() {
  projMatrix[0] = [0.601080953882, 0.335124580762, 0.0791121819909, 0.141895965096, 0.0224741130161, -0.499494658874, -0.243482260103, -0.436710875769];
  projMatrix[1] = [-0.0224741130161, 0.499494658874, 0.243482260103, 0.436710875769, 0.601080953882, 0.335124580762, 0.0791121819909, 0.141895965096];
  gramSchmidt();  // (Actually unnecessary)
}

// This class has order 18 and centralizer of order 108, and acts on the root by 13 cycles of length 18 and 3 cycles of length 2.
var rotate18 = [120, 121, 126, 122, 120, 132, 126, 122, 140, 150, 140, 132, 126, 122, 121, 162, 120, 122, 121, 126, 122, 120, 132, 126, 122, 121, 140, 132, 126, 122, 120, 150, 140, 132, 126, 122, 162, 120, 122, 121, 126, 122, 120, 132, 126, 122, 121, 162, 120, 122];

function chooseProjection18() {
  projMatrix[0] = [0.56082647979, 0.243265683415, 0.56082647979, 0.0359925942184, -0.0359925942184, -0.243265683415, -0.353553390593, -0.353553390593];
  projMatrix[1] = [0.0428943034666, 0.0988888398828, 0.365353986997, 0.668366972112, 0.260118681648, 0.507137130347, 0.204124145232, 0.204124145232];
  projMatrix[1][3] += 0.025;  // Blur!
  projMatrix[1][4] += 0.025;  // Blur!
  gramSchmidt();
}

// This class has order 14 and centralizer of order 28, and acts on the root by 17 cycles of length 14 and 1 cycle of length 2.
var rotate14 = [121, 122, 120, 121, 162, 120, 122, 121, 126, 122, 120, 132, 126, 122, 121, 140, 132, 126, 122, 120, 150, 140, 132, 162, 120, 122, 121, 126, 122, 120, 132, 126, 122, 121];

function chooseProjection14() {
  projMatrix[0] = [0.748849359032, 0.267261241912, 0.333269317529, 0.118942442321, -0.118942442321, -0.333269317529, -0.267261241912, -0.214326875208];
  projMatrix[1] = [0, 0.231920613924, 0.52112088917, 0.417906505941, 0.417906505941, 0.52112088917, 0.231920613924, 0];
  projMatrix[0][0] += 0.01;  // Blur!
  gramSchmidt();
}

function changeProjectionRandom() {
  for (var k = 0; k < 8; k++) {
    projMatrix[0][k] += gaussian() * 0.02;
    projMatrix[1][k] += gaussian() * 0.02;
  }
  gramSchmidt();
}

function customHyperPlane(x, y) {
  projMatrix[0] = x
  projMatrix[1] = y
  gramSchmidt(0)
}

function rotate(rotationMatrix) {
  for (var k = 0; k < 8; k++) {
    projMatrix[0][k] += rotationMatrix[k];
    projMatrix[1][k] += rotationMatrix[k];
  }
  gramSchmidt();
}

var rootProj = new Array();  // Projection of the roots (in canvas coordinates)
var rootDotSize = new Array();  // Dot sizes of the roots (in canvas pixels)
var rootColor = new Array();  // Root colors (RGB coordinates)

function computeProjections() {
  for (var n = 0; n < roots.length; n++) {
    // compute projections and resize by 100 and translated by 300
    rootProj[n] = [dotprod(projMatrix[0], roots[n]) * zoom + translateX,
    dotprod(projMatrix[1], roots[n]) * zoom + translateY];
  }
}

function computeColors() {
  for (var n = 0; n < roots.length; n++) {
    var x = dotprod(projMatrix[0], roots[n]) / 2.83;
    var y = dotprod(projMatrix[1], roots[n]) / 2.83;
    var hue = ((Math.atan2(-y, x) / Math.PI) + 1) * 3.;
    var red; var grn; var blu;
    if (hue < 1.) {
      red = 0.; grn = 1.; blu = hue;
    } else if (hue < 2.) {
      red = 0.; grn = 2. - hue; blu = 1.;
    } else if (hue < 3.) {
      red = hue - 2.; grn = 0.; blu = 1.;
    } else if (hue < 4.) {
      red = 1.; grn = 0.; blu = 4. - hue;
    } else if (hue < 5.) {
      red = 1.; grn = hue - 4.; blu = 0.;
    } else {
      red = 6. - hue; grn = 1.; blu = 0.;
    }
    var sat = Math.sqrt(x * x + y * y) * 0.9 + 0.1;
    red = sat * red + (1. - sat);
    grn = sat * grn + (1. - sat);
    blu = sat * blu + (1. - sat);
    rootColor[n] = [Math.floor(red * 255 + 0.5), Math.floor(grn * 255 + 0.5), Math.floor(blu * 255 + 0.5)];
  }
}

function computeDotSizes() {
  for (var n = 0; n < roots.length; n++) {
    var max = 0;
    for (var k = 0; k < 8; k++)
      if (Math.abs(roots[n][k]) > max)
        max = Math.abs(roots[n][k]);
    rootDotSize[n] = (max == 1) ? 5. : 7.;
  }
}

// The current Weyl group element, as a permutation of [0..239]
var permutation = new Array();

function initPerm() {
  for (var i = 0; i < roots.length; i++) {
    permutation[i] = i;
  }
}

function drawPoints() {  // Draw points using first color scheme.
  for (var n = 0; n < rootProj.length; n++) {
    var np = permutation[n];
    ctx1.fillStyle = "rgb(" + rootColor[np][0] + "," + rootColor[np][1] + "," + rootColor[np][2] + ")";
    ctx1.beginPath();
    ctx1.arc(
      rootProj[n][0], rootProj[n][1],
      rootDotSize[np],
      0, Math.PI * 2, true);
    ctx1.closePath();
    ctx1.fill();
    // ctx1.stroke()
  }
}


function clearCanvas() {
  ctx1.clearRect(0, 0, width, height);
  ctx0.clearRect(0, 0, width, height);
}


function drawLines() {  // Draw the 6720 edges
  ctx0.lineWidth = 0.1;
  ctx0.strokeStyle = "rgb(0,0,0)";
  for (var n = 0; n < rootProj.length; n++) {
    for (var m = 0; m < adjacent[n].length; m++) {
      var n2 = adjacent[n][m];
      if (n2 < n)
        continue;
      ctx0.beginPath();
      ctx0.moveTo(rootProj[n][0], rootProj[n][1]);
      ctx0.lineTo(rootProj[n2][0], rootProj[n2][1]);
      ctx0.stroke();
    }
  }
}

var intervalId
var speed = .03
function onLoad() {

  candiv = document.getElementById("candiv");
  canvas0 = document.getElementById("canvas0");
  canvas0.width = width
  canvas0.height = height
  canvas1 = document.getElementById("canvas1");
  canvas1.width = width
  canvas1.height = height
  if (typeof (canvas0.getContext) != "function") {
    alert("Your browser does not support the HTML5 <canvas> element.\n"
      + "This page will not function.");
    throw ("canvas unsupported");
  }

  ctx0 = canvas0.getContext("2d");
  ctx0.lineCap = "round"; ctx0.lineJoin = "round";
  ctx1 = canvas1.getContext("2d");
  ctx1.lineCap = "round"; ctx1.lineJoin = "round";
  createRoots();
  createAdjacent();
  chooseProjection30()
  paint()
  runInterval()
}

function paint() {
  clearCanvas()
  computeProjections()
  computeColors();
  computeDotSizes();
  initPerm();
  clearCanvas()
  drawLines();
  drawPoints();
}

var keys = {}
var playerspeed = 1
document.addEventListener("keydown", (e) => {
  keys[e.code] = "down"
})

document.addEventListener("keyup", (e) => {
  keys[e.code] = "up"
})


function rand(min, max) {
  return Math.random() * (max - min) + min;
}

var inputOn = false
function runInterval(rotateMatrix) {
  intervalId = setInterval(function () {
    if (!inputOn) {
      if (keys["ArrowLeft"] == "down") translateX += playerspeed * 10
      if (keys["ArrowRight"] == "down") translateX -= playerspeed * 10
      if (keys["ArrowUp"] == "down") translateY += playerspeed * 10
      if (keys["ArrowDown"] == "down") translateY -= playerspeed * 10
      if (keys["KeyQ"] == "down") zoom += playerspeed
      if (keys["KeyA"] == "down") zoom -= playerspeed

    }

    if (rotateMatrix) {
      const lastProj = JSON.parse(JSON.stringify(rootProj)).map(i => i.map(Math.floor))
      rotate(rotateMatrix)
      computeProjections();
      const newProj = JSON.parse(JSON.stringify(rootProj)).map(i => i.map(Math.floor))
      if (JSON.stringify(lastProj) == JSON.stringify(newProj)) {
        clearInterval(intervalId)
        runInterval([...Array(8).keys()].map(i => rand(-speed, speed)))
        return
      }
    }
    paint()
  }, 60)

}


class App extends Component {
  constructor(props) {
    super(props)
    this.x = JSON.stringify([2 - 4 / Math.sqrt(3), 0, 1 - 1 / Math.sqrt(3), 1 - 1 / Math.sqrt(3), 0, -1, 1, 0])
    this.y = JSON.stringify([0, -2 + 4 / Math.sqrt(3), -1 + 1 / Math.sqrt(3), 1 - 1 / Math.sqrt(3), 0, 1 / Math.sqrt(3), 1 / Math.sqrt(3), -2 / Math.sqrt(3)])

    this.state = {
      xtemp: "",
      ytemp: ""

    }
  }

  componentDidMount() {
    onLoad()
  }

  render() {

    return (
      <div className="App">
        <div style={{
          display: "flex",
          alignItems: "center",
          flexDirection: "column"
        }}>
          <div id="candiv" style={{ position: "relative", width: "100%", height: "100%" }}>
            <canvas id="canvas0" style={{ position: "absolute", left: 0 }} />
            <canvas id="canvas1" style={{ position: "absolute", left: 0 }} />
          </div>
        </div>
        <div style={{ position: "absolute", top: "0" }}>
          <pre className="mypre" style={{ padding: "1em", wordBreak: "break-word", whiteSpace: "pre-wrap" }}>
            E8, Gosset 4_21 polytope {"\n"}
            {roots.length} vertices {"\n"}
            8 dimensions {"\n\n"}

            made by Davide Andreazzini {"\n"}
            andreazzini.davide@gmail.com{"\n\n"}

            Use the Arrow Keys to navigate{"\n"}
            Use Q and A to zoom in and out{"\n\n"}
            Orthogonal plane coordinates {"\n"}
            x: {JSON.stringify(projMatrix.toJS()[0], "", "\t")}{"\n"}
            y: {JSON.stringify(projMatrix.toJS()[1], "", "\t")}
          </pre>
        </div>
        <div style={{ position: "absolute", right: 0 }}>
          <button onClick={() => {
            runInterval([...Array(8).keys()].map(i => rand(-speed, speed)))
          }}>start</button>
          <button onClick={() => {
            clearInterval(intervalId)
          }}>stop</button>
          <br />
          <button onClick={() => {
            chooseProjection30()
          }}>E8 coxeter plane (Petrie)</button>
          <br />

          <input type="text"
            onKeyDown={() => (inputOn = true)}
            onKeyUp={() => (inputOn = false)}

            onChange={e => {

              this.setState({ xtemp: e.target.value })

            }} /><br />
          <input type="text"
            onKeyDown={() => (inputOn = true)}
            onKeyUp={() => (inputOn = false)}

            onChange={e => {

              this.setState({ ytemp: e.target.value })
            }} /><br />
          <button onClick={() => {
            console.log(this.state)
            const x = JSON.parse(this.state.xtemp) || this.x
            const y = JSON.parse(this.state.ytemp) || this.y
            customHyperPlane(x, y)
          }}>custom hyperplane</button>

        </div>


      </div>
    );
  }
}

export default observer(App);
