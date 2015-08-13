var scene;
var aspect;
var camera;
var controls;
var renderer;
var worker;
var gui;

var plane;

var myTree;
var centers;
var tree;
var twig;


var light;
var lightDot;
var SunX=0;
var SunY=100;
var SunZ=0;

var contX;
var contY;
var contZ;

var preset='Top';
var maxGenerations=200;
var numOfInstances=12;
var numOfChild=8;
var sameBestIter=25;
var useBest=true;
var minFitnessTerm=0;
var mutationRate=0.2;

var trunkCollision=false;
var facing=false;

function workerMessage(evt){
	scene.remove(tree);
	scene.remove(twig);
	tree = treeGeometry(evt.data);
	twig = twigGeometry(evt.data);
	scene.add(tree);
	scene.add(twig);
}

this.Start=function(){
	var gaSettings=[maxGenerations,numOfInstances,numOfChild,sameBestIter,useBest,minFitnessTerm,mutationRate,facing,trunkCollision];
		//console.log("worker");
	worker.postMessage([1,gaSettings]);
}

function updateLight(){
	light.position.set( SunX, SunY, SunZ );//top
	var pos=[light.position.x,light.position.y,light.position.z];
	lightDot.position.set(light.position.x, light.position.y, light.position.z);
	worker.postMessage([0,pos]);
}

function init(){
	worker=new Worker("GAworker.js");
	worker.onmessage=workerMessage;
	scene = new THREE.Scene();
	scene.fog = new THREE.FogExp2( 0x9CCAF6, 0.05 );
	
	aspect = window.innerWidth / window.innerHeight;
	camera = new THREE.PerspectiveCamera( 75, aspect, 0.1, 1000 );
	renderer = new THREE.WebGLRenderer();
	renderer.setSize( window.innerWidth, window.innerHeight );
	renderer.shadowMapEnabled = true;
	renderer.shadowMapSoft = true;
	document.body.appendChild( renderer.domElement );
	camera.position.z = 10;
	camera.position.y = 3;
	
	controls = new THREE.OrbitControls( camera );
	controls.rotateSpeed = 1.0;
	controls.zoomSpeed = 1.2;
	controls.panSpeed = 0.8;
	controls.noZoom = false;
	controls.noPan = true;
	controls.staticMoving = true;
	controls.dynamicDampingFactor = 0.3;
	controls.target=new THREE.Vector3(0,3,0);

	renderer.shadowCameraNear = 3;
	renderer.shadowCameraFar = camera.far;
	renderer.shadowCameraFov = 50;

	renderer.shadowMapBias = 0.0039;
	renderer.shadowMapDarkness = 0.5;
	renderer.shadowMapWidth = 4096;
	renderer.shadowMapHeight = 4096;
	
	light = new THREE.PointLight( 0xffffff, 1 );
	//light.position.set( 50, 50, 50 );
	light.position.set( SunX, SunY, SunZ );//top
	//light.position.set( 100, 20, 0 );//side
	//light.position.set( 0, 20, 100 );
	scene.add( light );
	var pos=[light.position.x,light.position.y,light.position.z];
	worker.postMessage([0,pos]);
	
	var geo=new THREE.Geometry();
	geo.vertices.push(new THREE.Vector3( 0,  0, 0 ));
	pcMat = new THREE.PointCloudMaterial( { size: 15, sizeAttenuation: false, alphaTest: 0.5, transparent: true, fog: false} );
	lightDot=new THREE.PointCloud(geo,pcMat);
	lightDot.position.set(light.position.x, light.position.y, light.position.z);
	scene.add(lightDot);
	
	var directionalLight = new THREE.DirectionalLight( 0xffffff, 0.7 );
	directionalLight.position.set( 0, 1, 0 );
	scene.add( directionalLight );

	
	
	geometry = new THREE.BoxGeometry( 500, 0.01, 500 );
	planeMat = new THREE.MeshPhongMaterial( 
		{color: 0xEEEEEE, shading: THREE.SmoothShading}
		);
		 
	
	plane = new THREE.Mesh(geometry, planeMat);
	plane.receiveShadow = true;
	scene.add(plane);

	var sphere= new THREE.SphereGeometry(300, 32, 32);
	var skyMat = new THREE.MeshPhongMaterial( 
		{color: 0xFFFFFF, shading: THREE.SmoothShading}
		);
	skyMat.side = THREE.BackSide
	var sky=new THREE.Mesh(sphere,skyMat);
	scene.add(sky);
	
	myTree=randomTree();
	
	tree = treeGeometry(myTree);
	scene.add( tree );
	twig = twigGeometry(myTree);
	/*
	centers=leafCenters(myTree);
	normals=leafNormals(myTree);
	var geo=new THREE.Geometry(); //leafCenters,normal debug
	for(i=0;i<centers.length;i++){//centers.length
			geo.vertices.push(new THREE.Vector3( centers[i][0],  centers[i][1], centers[i][2] ));
			geo.vertices.push(new THREE.Vector3( centers[i][0]+normals[i][0]/3,  centers[i][1]+normals[i][1]/3, centers[i][2]+normals[i][2]/3 ));
	}
	*/
	/*var geo=treeLightFitness(myTree);
	pcMat = new THREE.PointCloudMaterial( { size: 5, sizeAttenuation: false, alphaTest: 0.5, transparent: true } );
	twigNormals=new THREE.PointCloud(geo,pcMat);
	scene.add(twigNormals);*/
	
	scene.add(twig);
	gui = new dat.GUI();
	var f1 = gui.addFolder("Sun position");
	//var cont=f1.add(self, 'preset', [ 'Top', 'Side1', 'Side2' ] );
	
	var contX=f1.add(self, 'SunX').step(1);	
	contX.onChange(updateLight);
	
	var contY=f1.add(self, 'SunY').step(1);	
	contY.onChange(updateLight);
	
	var contZ=f1.add(self, 'SunZ').step(1);	
	contZ.onChange(updateLight);
	
	var f2 = gui.addFolder("Genetic Algorithm");
	var f3=f2.addFolder("Fitness settings");
	f2.add(self, 'maxGenerations').min(0).step(1);
	f2.add(self, 'numOfInstances').min(0).step(1);
	f2.add(self, 'numOfChild').min(0).step(2);
	f2.add(self, 'sameBestIter').min(0).step(1);
	f2.add(self, 'useBest');
	f2.add(self, 'minFitnessTerm').min(0).step(1);
	f2.add(self,'mutationRate').min(0).max(1).step(0.1);
	f2.add(self, 'Start');
	
	f3.add(self,'trunkCollision');
	f3.add(self,'facing');
	
	
}


var render = function () {
	requestAnimationFrame( render );
	
	if(keyboard.pressed("r") == true)	{
		myTree=randomTree();	
		
		scene.remove(tree);
		scene.remove(twig);
		tree = treeGeometry(myTree);
		twig = twigGeometry(myTree);
		scene.add(tree);
		scene.add(twig);
		
	}
	
	controls.update();
	renderer.render( scene, camera );
};

window.addEventListener("resize", function()
{
	renderer.setSize(window.innerWidth, window.innerHeight);
	camera.aspect = window.innerWidth / window.innerHeight;
	camera.updateProjectionMatrix();
});

function randomTree(){
	var seed=Math.round(Math.random()*4000);
	var initalBranchLength=0.5+Math.random()*0.3;
	var lengthFalloffFactor=0.5+Math.random()*0.3;
	var lengthFalloffPower=0.3+Math.random()*0.4;
	var clumpMax=0.4+Math.random()*0.1;
	var clumpMin=clumpMax-Math.random()*0.4;
	var branchFactor=2.0+Math.random()*2.0;
	var dropAmount=-0.3+Math.random()*0.6;
	var growAmount=-0.5+Math.random()*1.5;
	var sweepAmount=-0.05+Math.random()*0.1;
	var climbRate=0.05+Math.random()*0.95;
	var trunkKink=Math.random()*0.3;
	var taperRate=0.7+Math.random()*0.3;
	var radiusFalloffRate=0.74+Math.random()*0.05;
	var twistRate=Math.random()*10.0;
	var trunkLength=1.5+Math.random()*1.6;
	var myTree = new Tree(
		{"seed":seed,
		"segments":6,
		"levels":5,
		"vMultiplier":1.16,
		"twigScale":0.2,
		"initalBranchLength":initalBranchLength,
		"lengthFalloffFactor":lengthFalloffFactor,
		"lengthFalloffPower":lengthFalloffPower,
		"clumpMax":clumpMax,
		"clumpMin":clumpMin,
		"branchFactor":branchFactor,
		"dropAmount":dropAmount,
		"growAmount":growAmount,
		"sweepAmount":sweepAmount,
		"maxRadius":0.111,
		"climbRate":climbRate,
		"trunkKink":trunkKink,
		"treeSteps":4,
		"taperRate":taperRate,
		"radiusFalloffRate":radiusFalloffRate,
		"twistRate":twistRate,
		"trunkLength":trunkLength,
		"trunkMaterial":"TrunkType1",
		"twigMaterial":"BranchType6"}
		
	);
	
	myTree.seed=seed;
	myTree.initalBranchLength=initalBranchLength;
	myTree.lengthFalloffFactor=lengthFalloffFactor;
	myTree.lengthFalloffPower=lengthFalloffPower;
	myTree.clumpMax=clumpMax;
	myTree.clumpMin=clumpMin;
	myTree.branchFactor=branchFactor;
	myTree.dropAmount=dropAmount;
	myTree.growAmount=growAmount;
	myTree.sweepAmount=sweepAmount;
	myTree.climbRate=climbRate;
	myTree.trunkKink=trunkKink;
	myTree.taperRate=taperRate;
	myTree.radiusFalloffRate=radiusFalloffRate;
	myTree.twistRate=twistRate;
	myTree.trunkLength=trunkLength;
	//console.log("verts: "+myTree.verts.length+", polygons: "+myTree.faces.length+", twigPolygons: "+myTree.facesTwig.length);
	
	return myTree;
}

function treeGeometry(myTree){
	var treeGeo=new THREE.Geometry();
	
	var treeMat = new THREE.MeshPhongMaterial( 
	{color: 0xCCCCCC,//0x804515,
	 shading: THREE.FlatShading
	 } 
	);
	for(i=0;i<myTree.verts.length;i++){
		treeGeo.vertices.push(new THREE.Vector3( myTree.verts[i][0],  myTree.verts[i][1], myTree.verts[i][2] ));
	}
	
	for(i=0;i<myTree.faces.length;i++){
		treeGeo.faces.push(new THREE.Face3( myTree.faces[i][0],  myTree.faces[i][1], myTree.faces[i][2] ));
	}
	
	treeGeo.computeBoundingSphere();
	treeGeo.computeFaceNormals ();
	treeGeo.computeVertexNormals ();
	var tree = new THREE.Mesh( treeGeo, treeMat );
	tree.castShadow = true;
	tree.receiveShadow = true;
	return tree;
}

function twigGeometry(myTree){
	var twigGeo=new THREE.Geometry();
	var twigMat = new THREE.MeshPhongMaterial( 
	{color: 0xFFFFFF,//0x84E533,
	 shading: THREE.SmoothShading} 
	);
	for(i=0;i<myTree.vertsTwig.length;i++){//myTree.vertsTwig.length
		twigGeo.vertices.push(new THREE.Vector3( myTree.vertsTwig[i][0],  myTree.vertsTwig[i][1], myTree.vertsTwig[i][2] ));
	}
	
	for(i=0;i<myTree.facesTwig.length;i++){//myTree.facesTwig.length
		twigGeo.faces.push(new THREE.Face3( myTree.facesTwig[i][0],  myTree.facesTwig[i][1], myTree.facesTwig[i][2] ));
	}
	
	twigGeo.computeBoundingSphere();
	twigGeo.computeFaceNormals ();
	var twig = new THREE.Mesh( twigGeo, twigMat );
	twig.castShadow = true;
	twig.receiveShadow = true;
	return twig;
}
	
var keyboard = new THREEx.KeyboardState();


init();
render();