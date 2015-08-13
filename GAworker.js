var window=self;

importScripts('threejs/build/three.min.js'); 
importScripts('lib/proctree.js'); 
var light = new THREE.PointLight( 0xffffff, 1 );
//light.position.set( 50, 50, 50 );
//light.position.set( 0, 500, 0 );//top
//light.position.set( 100, 20, 0 );//side
light.position.set( 0, 20, 100 );
var working=false;



onmessage=serviceMessage;

function serviceMessage(evt){
	//var sReturnMessage = ("Hello from the worker thread! This is what you sent my way: " + evt.data); 
	//console.log(evt.data[0]);
	//console.log(evt.data[1]);
	
	var msg=evt.data[0];
	var args=evt.data[1];
	if(msg==1){
		if(working){
			console.log("Worker is busy");
		}else{
			var a=GAOptimize(args[0],args[1],args[2],args[3],args[4],args[5],args[6],args[7],args[8]);
			//console.log(args);
			postMessage(a);
		}
	}else if(msg==0){
		
		light.position.set(args[0],args[1],args[2]);
		//console.log("Position set: "+light.position.x+", "+light.position.y+", "+light.position.z);
	}
	// Pass our message back to the creator's thread 
	
}


var dot=function(v1,v2){
	return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
};

var cross=function(v1,v2){
	return [v1[1]*v2[2]-v1[2]*v2[1],v1[2]*v2[0]-v1[0]*v2[2],v1[0]*v2[1]-v1[1]*v2[0]];
};

var length=function(v){
	return Math.sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
};

var normalize=function(v){
	var l=length(v);
	return scaleVec(v,1/l);
};

var scaleVec=function(v,s){
	return [v[0]*s,v[1]*s,v[2]*s];
};

var subVec=function(v1,v2){
	return [v1[0]-v2[0],v1[1]-v2[1],v1[2]-v2[2]];
};

var addVec=function(v1,v2){
	return [v1[0]+v2[0],v1[1]+v2[1],v1[2]+v2[2]];
};

var distance=function(v1,v2){
	return Math.sqrt(Math.pow((v1[0]-v2[0]),2)+Math.pow((v1[1]-v2[1]),2)+Math.pow((v1[2]-v2[2]),2));
};

var sum=function(v1){
	var s=0;
	
	for(i=0;i<v1.length;i++){
		s=s+v1[i];
	}
	return s;
};

function leafNormals(myTree){
	var normals=new Array(myTree.facesTwig.length);
	for(i=0;i<myTree.facesTwig.length;i++){
		var face=myTree.facesTwig[i];
		//console.log(face);
		var u=subVec(myTree.vertsTwig[face[1]],myTree.vertsTwig[face[0]]);
		var v=subVec(myTree.vertsTwig[face[2]],myTree.vertsTwig[face[0]]);
		var norm=normalize(cross(u,v));
		normals[i]=[norm[0],norm[1],norm[2]];
	}
	return normals;
}

function trunkNormals(myTree){
	var normals=new Array(myTree.faces.length);
	for(i=0;i<myTree.faces.length;i++){
		var face=myTree.faces[i];
		//console.log(face);
		var u=subVec(myTree.verts[face[1]],myTree.verts[face[0]]);
		var v=subVec(myTree.verts[face[2]],myTree.verts[face[0]]);
		var norm=normalize(cross(u,v));
		normals[i]=[norm[0],norm[1],norm[2]];
	}
	return normals;
}


function treeLightFitness(myTree){
	var fit=0;
	var p0=[light.position.x,light.position.y,light.position.z]
	var centers=leafCenters(myTree);
	var normals=leafNormals(myTree);
	var fits=new Array(normals.length);
	
	var geo=new THREE.Geometry();
	
	for(i=0; i<normals.length; i++){
		fits[i]=1;
		var p=centers[i];
		var n=normals[i];
		
		var d=-1*dot(p,n);
		var v=normalize(subVec(p,p0));
		var t0=-1*(dot(p0,n)+d)/dot(v,n);//dela pravilno
		var P=addVec(p0,scaleVec(v,t0));
		//console.log(p);
		//console.log(P);
		//console.log(t0);
		for(j=0;j<normals.length;j++){
			if(j==i){
				continue;
			}
			
			var dist=distance(centers[i],centers[j]);
			
			if(dist<0.001){
				continue;
			}
			
			p=centers[j];//plane poly2
			n=normals[j];
			d=-1*dot(p,n);
			
			var t=-1*(dot(p0,n)+d)/dot(v,n);
			P=addVec(p0,scaleVec(v,t));
			
			//console.log(Math.round(dot(P,n)+d));//dela pravilno
			if(t<t0){
				var dif=Math.abs(t-t0);
				if(dif>0.001){
					
					var A=myTree.vertsTwig[myTree.facesTwig[j][0]];
					var B=myTree.vertsTwig[myTree.facesTwig[j][1]];
					var C=myTree.vertsTwig[myTree.facesTwig[j][2]];
					//inTriangle3(A,B,C,P);
					//console.log(A);
					if(nearTriangle(A,p,P)){
						if(inTriangle3(A,B,C,P)){//inTriangle2(A,B,C,P)
							fits[i]=0;
							/*geo.vertices.push(new THREE.Vector3( A[0],  A[1], A[2]));
							geo.vertices.push(new THREE.Vector3( B[0],  B[1], B[2]));
							geo.vertices.push(new THREE.Vector3( C[0],  C[1], C[2]));*/
							//geo.vertices.push(new THREE.Vector3( P[0],  P[1], P[2]));
							//console.log(p+", "+P);
							break;
						}
					}
				}
			}
		}
	}

	for(i=0;i<fits.length;i++){
		fit=fit+fits[i];
	}
	//return geo;
	return fit;
}

function treeLightFacingFitness(myTree){

	var fit=0;
	var p0=[light.position.x,light.position.y,light.position.z]
	var centers=leafCenters(myTree);
	var normals=leafNormals(myTree);
	/*if(isNaN(normals[0])){
		console.log(normals[0]);
	}*/
	//console.log(normals.length);
	var fits=new Array(normals.length);
	
	//var geo=new THREE.Geometry();
	
	for(i=0; i<normals.length; i++){
		fits[i]=1;
		var p=centers[i];
		var n=normals[i];
		//console.log(p);
		
		var d=-1*dot(p,n);
		var v=normalize(subVec(p,p0));
		var t0=-1*(dot(p0,n)+d)/dot(v,n);//dela pravilno
		var P=addVec(p0,scaleVec(v,t0));
		//console.log(p);
		//console.log(P);
		//console.log(t0);
		for(j=0;j<normals.length;j++){
			if(j==i){
				continue;
			}
			
			var dist=distance(centers[i],centers[j]);
			
			if(dist<0.001){
				continue;
			}
			
			p=centers[j];//plane poly2
			n=normals[j];
			
			fits[i]=Math.abs(dot(v,n));
			
			d=-1*dot(p,n);
			
			var t=-1*(dot(p0,n)+d)/dot(v,n);
			P=addVec(p0,scaleVec(v,t));
			
			//console.log(Math.round(dot(P,n)+d));//dela pravilno
			if(t<t0){
				var dif=Math.abs(t-t0);
				if(dif>0.001){
					
					var A=myTree.vertsTwig[myTree.facesTwig[j][0]];
					var B=myTree.vertsTwig[myTree.facesTwig[j][1]];
					var C=myTree.vertsTwig[myTree.facesTwig[j][2]];
					//inTriangle3(A,B,C,P);
					//console.log(A);
					if(nearTriangle(A,p,P)){
						if(inTriangle3(A,B,C,P)){//inTriangle2(A,B,C,P)
							fits[i]=0;
							/*geo.vertices.push(new THREE.Vector3( A[0],  A[1], A[2]));
							geo.vertices.push(new THREE.Vector3( B[0],  B[1], B[2]));
							geo.vertices.push(new THREE.Vector3( C[0],  C[1], C[2]));
							geo.vertices.push(new THREE.Vector3( P[0],  P[1], P[2]));*/
							//console.log(p+", "+P);
							break;
						}
					}
				}
			}
		}
	}

	for(i=0;i<fits.length;i++){
		fit=fit+fits[i];
		/*if(isNaN(fits[i])){
			console.log(myTree);
		}*/
	}
	//return geo;
	return fit;
}

function treeTrunkLightFacingFitness(myTree){

	var fit=0;
	var p0=[light.position.x,light.position.y,light.position.z]
	var centers=leafCenters(myTree);
	var normals=leafNormals(myTree);
	var tcenters=trunkCenters(myTree);
	var tnormals=trunkNormals(myTree);
	/*if(isNaN(normals[0])){
		console.log(normals[0]);
	}*/
	//console.log(normals.length);
	var fits=new Array(normals.length);
	
	//var geo=new THREE.Geometry();
	
	for(i=0; i<normals.length; i++){
		fits[i]=1;
		var p=centers[i];
		var n=normals[i];
		//console.log(p);
		
		var d=-1*dot(p,n);
		var v=normalize(subVec(p,p0));
		var t0=-1*(dot(p0,n)+d)/dot(v,n);//dela pravilno
		var P=addVec(p0,scaleVec(v,t0));
		//console.log(p);
		//console.log(P);
		//console.log(t0);
		for(j=0;j<normals.length;j++){
			if(j==i){
				continue;
			}
			
			var dist=distance(centers[i],centers[j]);
			
			if(dist<0.001){
				continue;
			}
			
			p=centers[j];//plane poly2
			n=normals[j];
			
			fits[i]=Math.abs(dot(v,n));
			
			d=-1*dot(p,n);
			
			var t=-1*(dot(p0,n)+d)/dot(v,n);
			P=addVec(p0,scaleVec(v,t));
			
			//console.log(Math.round(dot(P,n)+d));//dela pravilno
			if(t<t0){
				var dif=Math.abs(t-t0);
				if(dif>0.001){
					
					var A=myTree.vertsTwig[myTree.facesTwig[j][0]];
					var B=myTree.vertsTwig[myTree.facesTwig[j][1]];
					var C=myTree.vertsTwig[myTree.facesTwig[j][2]];
					//inTriangle3(A,B,C,P);
					//console.log(A);
					if(nearTriangle(A,p,P)){
						if(inTriangle3(A,B,C,P)){//inTriangle2(A,B,C,P)
							fits[i]=0;
							/*geo.vertices.push(new THREE.Vector3( A[0],  A[1], A[2]));
							geo.vertices.push(new THREE.Vector3( B[0],  B[1], B[2]));
							geo.vertices.push(new THREE.Vector3( C[0],  C[1], C[2]));
							geo.vertices.push(new THREE.Vector3( P[0],  P[1], P[2]));*/
							//console.log(p+", "+P);
							break;
						}
					}
				}
			}
		}
		if(fits[i]!=0){
			for(j=0;j<tnormals.length;j++){
				p=tcenters[j];//plane poly2
				n=tnormals[j];
				
				fits[i]=Math.abs(dot(v,n));
				
				d=-1*dot(p,n);
				
				var t=-1*(dot(p0,n)+d)/dot(v,n);
				P=addVec(p0,scaleVec(v,t));
				
				//console.log(Math.round(dot(P,n)+d));//dela pravilno
				if(t<t0){
					var dif=Math.abs(t-t0);
					if(dif>0.001){
						
						var A=myTree.verts[myTree.faces[j][0]];
						var B=myTree.verts[myTree.faces[j][1]];
						var C=myTree.verts[myTree.faces[j][2]];
						//inTriangle3(A,B,C,P);
						//console.log(A);
						if(nearTriangle(A,p,P)){
							if(inTriangle3(A,B,C,P)){//inTriangle2(A,B,C,P)
								fits[i]=0;
								/*geo.vertices.push(new THREE.Vector3( A[0],  A[1], A[2]));
								geo.vertices.push(new THREE.Vector3( B[0],  B[1], B[2]));
								geo.vertices.push(new THREE.Vector3( C[0],  C[1], C[2]));
								geo.vertices.push(new THREE.Vector3( P[0],  P[1], P[2]));*/
								//console.log(p+", "+P);
								break;
							}
						}
					}
				}
			}
		}
	}

	for(i=0;i<fits.length;i++){
		fit=fit+fits[i];
		/*if(isNaN(fits[i])){
			console.log(myTree);
		}*/
	}
	//return geo;
	return fit;
}

function treeTrunkLightFitness(myTree){

	var fit=0;
	var p0=[light.position.x,light.position.y,light.position.z]
	var centers=leafCenters(myTree);
	var normals=leafNormals(myTree);
	var tcenters=trunkCenters(myTree);
	var tnormals=trunkNormals(myTree);
	/*if(isNaN(normals[0])){
		console.log(normals[0]);
	}*/
	//console.log(normals.length);
	var fits=new Array(normals.length);
	
	//var geo=new THREE.Geometry();
	
	for(i=0; i<normals.length; i++){
		fits[i]=1;
		var p=centers[i];
		var n=normals[i];
		//console.log(p);
		
		var d=-1*dot(p,n);
		var v=normalize(subVec(p,p0));
		var t0=-1*(dot(p0,n)+d)/dot(v,n);//dela pravilno
		var P=addVec(p0,scaleVec(v,t0));
		//console.log(p);
		//console.log(P);
		//console.log(t0);
		for(j=0;j<normals.length;j++){
			if(j==i){
				continue;
			}
			
			var dist=distance(centers[i],centers[j]);
			
			if(dist<0.001){
				continue;
			}
			
			p=centers[j];//plane poly2
			n=normals[j];
			
			//fits[i]=Math.abs(dot(v,n));
			
			d=-1*dot(p,n);
			
			var t=-1*(dot(p0,n)+d)/dot(v,n);
			P=addVec(p0,scaleVec(v,t));
			
			//console.log(Math.round(dot(P,n)+d));//dela pravilno
			if(t<t0){
				var dif=Math.abs(t-t0);
				if(dif>0.001){
					
					var A=myTree.vertsTwig[myTree.facesTwig[j][0]];
					var B=myTree.vertsTwig[myTree.facesTwig[j][1]];
					var C=myTree.vertsTwig[myTree.facesTwig[j][2]];
					//inTriangle3(A,B,C,P);
					//console.log(A);
					if(nearTriangle(A,p,P)){
						if(inTriangle3(A,B,C,P)){//inTriangle2(A,B,C,P)
							fits[i]=0;
							/*geo.vertices.push(new THREE.Vector3( A[0],  A[1], A[2]));
							geo.vertices.push(new THREE.Vector3( B[0],  B[1], B[2]));
							geo.vertices.push(new THREE.Vector3( C[0],  C[1], C[2]));
							geo.vertices.push(new THREE.Vector3( P[0],  P[1], P[2]));*/
							//console.log(p+", "+P);
							break;
						}
					}
				}
			}
		}
		if(fits[i]!=0){
			for(j=0;j<tnormals.length;j++){
				p=tcenters[j];//plane poly2
				n=tnormals[j];
				
				//fits[i]=Math.abs(dot(v,n));
				
				d=-1*dot(p,n);
				
				var t=-1*(dot(p0,n)+d)/dot(v,n);
				P=addVec(p0,scaleVec(v,t));
				
				//console.log(Math.round(dot(P,n)+d));//dela pravilno
				if(t<t0){
					var dif=Math.abs(t-t0);
					if(dif>0.001){
						
						var A=myTree.verts[myTree.faces[j][0]];
						var B=myTree.verts[myTree.faces[j][1]];
						var C=myTree.verts[myTree.faces[j][2]];
						//inTriangle3(A,B,C,P);
						//console.log(A);
						if(nearTriangle(A,p,P)){
							if(inTriangle3(A,B,C,P)){//inTriangle2(A,B,C,P)
								fits[i]=0;
								/*geo.vertices.push(new THREE.Vector3( A[0],  A[1], A[2]));
								geo.vertices.push(new THREE.Vector3( B[0],  B[1], B[2]));
								geo.vertices.push(new THREE.Vector3( C[0],  C[1], C[2]));
								geo.vertices.push(new THREE.Vector3( P[0],  P[1], P[2]));*/
								//console.log(p+", "+P);
								break;
							}
						}
					}
				}
			}
		}
	}

	for(i=0;i<fits.length;i++){
		fit=fit+fits[i];
		/*if(isNaN(fits[i])){
			console.log(myTree);
		}*/
	}
	//return geo;
	return fit;
}

function inTriangle(A,B,C,P){
	areaABC=length(cross(subVec(B,A),subVec(C,A)))/2;
	
	alpha=length(cross(subVec(B,P),subVec(C,P)))/(2*areaABC);
	beta=length(cross(subVec(C,P),subVec(A,P)))/(2*areaABC);
	gama=1-alpha-beta;
	
	if(alpha<0 || alpha>1){
		return false;
	}else if(beta<0 || beta>1){
		return false;
	}else if(gama<0 || gama>1){
		return false;
	}
	return true;
}

function inTriangle2(A,B,C,P){
	var v0=subVec(C,A);
	var v1=subVec(B,A);
	var v2=subVec(P,A);
	
	var dot00=dot(v0,v0);
	var dot01=dot(v0,v2);
	var dot02=dot(v0,v2);
	var dot11=dot(v1,v1);
	var dot12=dot(v1,v2);
	
	var invDenom = 1 /(dot00 * dot11 - dot01*dot01);
	var u = (dot11 * dot02 - dot01 * dot12) * invDenom;
	var v = (dot00 * dot12 - dot01 * dot02) * invDenom;

	// Check if point is in triangle
	return (u >= 0) && (v >= 0) && (u + v < 1)

}

function inTriangle3(A,B,C,P){
	var v0=subVec(A,P);
	var v1=subVec(B,P);
	var v2=subVec(C,P);

	var angAB=Math.acos((dot(v0,v1)/(length(v0)*length(v1))));
	var angBC=Math.acos((dot(v1,v2)/(length(v1)*length(v2))));
	var angCA=Math.acos((dot(v2,v0)/(length(v2)*length(v0))));
	
	var comb=angAB+angBC+angCA;
	
	if(comb<(Math.PI*2-0.01)){
		return false;
	}else{
		return true;
	}
}

function nearTriangle(A,center,P){
	var maxleng=distance(A,center);
	var pCentDist=distance(P,center);
	//console.log(maxleng);
	//console.log(pCentDist);
	if(pCentDist>maxleng){
		return false;
	}
	return true;
}

function leafCenters(myTree){
	var centers=new Array(myTree.facesTwig.length);
	for(i=0;i<myTree.facesTwig.length;i++){
		
		centers[i]=[(myTree.vertsTwig[ myTree.facesTwig[i][0] ][0]+myTree.vertsTwig[ myTree.facesTwig[i][1] ][0]+myTree.vertsTwig[ myTree.facesTwig[i][2] ][0])/3,
						 (myTree.vertsTwig[ myTree.facesTwig[i][0] ][1]+myTree.vertsTwig[ myTree.facesTwig[i][1] ][1]+myTree.vertsTwig[ myTree.facesTwig[i][2] ][1])/3,
						 (myTree.vertsTwig[ myTree.facesTwig[i][0] ][2]+myTree.vertsTwig[ myTree.facesTwig[i][1] ][2]+myTree.vertsTwig[ myTree.facesTwig[i][2] ][2])/3];
		/*if(isNaN(centers[i])){
			console.log(centers[i]);
		}*/
	}
	return centers;
}

function trunkCenters(myTree){
	var centers=new Array(myTree.faces.length);
	for(i=0;i<myTree.faces.length;i++){
		
		centers[i]=[(myTree.verts[ myTree.faces[i][0] ][0]+myTree.verts[ myTree.faces[i][1] ][0]+myTree.verts[ myTree.faces[i][2] ][0])/3,
						 (myTree.verts[ myTree.faces[i][0] ][1]+myTree.verts[ myTree.faces[i][1] ][1]+myTree.verts[ myTree.faces[i][2] ][1])/3,
						 (myTree.verts[ myTree.faces[i][0] ][2]+myTree.verts[ myTree.faces[i][1] ][2]+myTree.verts[ myTree.faces[i][2] ][2])/3];
		/*if(isNaN(centers[i])){
			console.log(centers[i]);
		}*/
	}
	return centers;
}

function GAOptimize(maxIter,numOfInstances,numOfChild,sameBestIter,useBest,minFitnessTerm, mutationRate, useFacing, useTrunk){
	working=true;

	if (typeof(minFitnessTerm)==='undefined'){
		minFitnessTerm=0;
	}
	
	console.log("Starting genetic algorithem: maxIter: "+maxIter+" numOfChild: "+numOfChild+" useBest: "+useBest+" minFitnessTerm: "+minFitnessTerm);
	
	var bestTree=randomTree();
	bestTree.fitness=0;
	
	/*var secondBest=randomTree();
	secondBest.fitness=0;*/
	
	var trees=[];
	var fitness=[];
	var newTrees=[];
	var bestIter=0;
	
	for( var i=0;i<numOfInstances;i++){//init
		var tree=randomTree();
		//treeLightFacingFitness(tree);
		//console.log(fit);
		trees.push(tree);
		
	}
	
	for(iter=0;iter<maxIter;iter++){
		if(useBest){
			trees.push(bestTree);
		}
		if(sameBestIter<=bestIter){
			console.log("Same best for "+sameBestIter+" iterations");
			break;
		}
		
		if(bestTree.fitness>minFitnessTerm && minFitnessTerm!=0){
			console.log("Terimantion fitness reached");
			break;
		}
		//console.log(trees);
		
		for(i=0;i<trees.length;i++){
			var fit=0;
			if(useFacing){
				if(useTrunk){
					fit=treeTrunkLightFacingFitness(trees[i]);
				}else{
					fit=treeLightFacingFitness(trees[i]);
				}
			}else{
				if(useTrunk){
					fit=treeTrunkLightFitness(trees[i]);
				}else{
					fit=treeLightFitness(trees[i]);
				}
			}
			//var fit=treeTrunkLightFitness(trees[i])//treeTrunkLightFacingFitness(trees[i]);
			if(isNaN(fit)){
				console.log(fit);
				fit=0;
				
			}
			fitness[i]=fit;
		}
		
		
		var all=[];//sort
		
		for(var i=0; i<trees.length;i++){
			all.push([ fitness[i], trees[i] ]);
		}
		
		all.sort(function(a, b) {
			return b[0] - a[0];
		});
		trees=[];
		fitness=[];
		
		for (var i = 0; i < all.length; i++) {
		   fitness.push(all[i][0]);
		   trees.push(all[i][1]);
		} 
		
		console.log(fitness);
		
		bestIter++;
		
		if(bestTree.fitness<fitness[0]){
		
			bestTree=trees[0];
			bestTree.fitness=fitness[0];
			bestIter=0;
			
			postMessage(bestTree);
		}
		
		console.log(iter+": Best Tree: "+bestTree.fitness);
		
		var selection=stohasticSelection(fitness,numOfChild);
		//console.log(selection);
		while(selection.length>0){
			console.log(selection);
			var sel1=selection[0];
			var b=1;
			var sel2=selection[b];
			
			if(selection.length==2){
				sel1=selection[0];
				sel2=selection[1];
				b=1;
			}else{
			
				while(fitness[sel1]==fitness[sel2]){
					b++;
					/*if(b>=selection.length){
						b--;
						break;
					}*/
					sel2=selection[b];
					if (typeof(sel2)==='undefined'){
						b=b-1;
						sel2=selection[b];
						break;
					}
				}
				
				while(sel1==sel2){
					b++;
					/*if(b>=selection.length){
						b--;
						break;
					}*/
					sel2=selection[b];
					if (typeof(sel2)==='undefined'){
						b=b-1;
						sel2=selection[b];
						break;
					}
				}

			}
			console.log(sel1+" "+sel2);
			var tempTrees=crossOver(trees[sel1],trees[sel2]);
			
			//console.log(tempTrees);
			//console.log(randomTree());
			for(i=0;i<tempTrees.length;i++){
				var mutationRand=Math.random();
				
				if(mutationRand<mutationRate){
					tempTrees[i]=mutateTree(tempTrees[i], 0.2)
				}
				
				newTrees.push(tempTrees[i]);
			}
			
			selection.splice(0,1);
			selection.splice(b-1,1);
		}
		
		for(i=newTrees.length;i<numOfInstances;i++){
			var tree=randomTree();
			newTrees.push(tree);
		}
		
		trees=[];
		for(i=0;i<newTrees.length;i++){
			trees[i]=newTrees[i];
		}
		//console.log(trees);
		newTrees=[];
		//console.log(newTrees.length);
		//var newTree=crossOver(trees[0],trees[1]);
		//console.log(trees[0]);
	}
	working=false;
	return bestTree;
	
}

function stohasticSelection(fitness,n){
	//console.log(fitness);
	var f=sum(fitness);
	var p=f/n;
	var curr=Math.random()*p;
	//console.log(f);
	var pointers=[];
	for(i=0;i<n;i++){
		var s=0;
		for(j=0;j<fitness.length;j++){
			s=s+fitness[j];
			if(curr<s){
				pointers[i]=j;
				
				break;
			}
			
		}
		curr+=p;
	}
	return pointers;
}

function crossOver(tree1,tree2){//enaka drevesa?
	var seed=[tree1.seed, tree2.seed];
	//console.log(seed);
	var initalBranchLength=[tree1.initalBranchLength, tree2.initalBranchLength];
	//console.log(initalBranchLength);
	var lengthFalloffFactor=[tree1.lengthFalloffFactor, tree2.lengthFalloffFactor];
	//console.log(lengthFalloffFactor);
	var lengthFalloffPower=[tree1.lengthFalloffPower, tree2.lengthFalloffPower];
	//console.log(lengthFalloffPower);
	var clumpMax=[tree1.clumpMax, tree2.clumpMax];
	//console.log(clumpMax);
	var clumpMin=[tree1.clumpMin, tree2.clumpMin];
	//console.log(clumpMin);
	var branchFactor=[tree1.branchFactor, tree2.branchFactor];
	//console.log(branchFactor);
	var dropAmount=[tree1.dropAmount, tree2.dropAmount];
	//console.log(dropAmount);
	var growAmount=[tree1.growAmount, tree2.growAmount];
	//console.log(growAmount);
	var sweepAmount=[tree1.sweepAmount, tree2.sweepAmount];
	//console.log(sweepAmount);
	var climbRate=[tree1.climbRate, tree2.climbRate];
	//console.log(climbRate);
	var trunkKink=[tree1.trunkKink, tree2.trunkKink];
	//console.log(trunkKink);
	var taperRate=[tree1.taperRate, tree2.taperRate];
	//console.log(taperRate);
	var radiusFalloffRate=[tree1.radiusFalloffRate, tree2.radiusFalloffRate];
	//console.log(radiusFalloffRate);
	var twistRate=[tree1.twistRate, tree2.twistRate];
	//console.log(twistRate);
	var trunkLength=[tree1.trunkLength, tree2.trunkLength];
	//console.log(trunkLength);
	
	
	var crossSel=[];
	for (i=0;i<16;i++){
		var ran=Math.random();
		if(ran>0.5){
			crossSel[i]=1;
		}else{
			crossSel[i]=0;
		}
	}
	//console.log(crossSel);
	var newTree1 = new Tree(
		{"seed":seed[crossSel[0]],
		"segments":6,
		"levels":5,
		"vMultiplier":1.16,
		"twigScale":0.2,
		"initalBranchLength":initalBranchLength[crossSel[1]],
		"lengthFalloffFactor":lengthFalloffFactor[crossSel[2]],
		"lengthFalloffPower":lengthFalloffPower[crossSel[3]],
		"clumpMax":clumpMax[crossSel[4]],
		"clumpMin":clumpMin[crossSel[5]],
		"branchFactor":branchFactor[crossSel[6]],
		"dropAmount":dropAmount[crossSel[7]],
		"growAmount":growAmount[crossSel[8]],
		"sweepAmount":sweepAmount[crossSel[9]],
		"maxRadius":0.111,
		"climbRate":climbRate[crossSel[10]],
		"trunkKink":trunkKink[crossSel[11]],
		"treeSteps":4,
		"taperRate":taperRate[crossSel[12]],
		"radiusFalloffRate":radiusFalloffRate[crossSel[13]],
		"twistRate":twistRate[crossSel[14]],
		"trunkLength":trunkLength[crossSel[15]],
		"trunkMaterial":"TrunkType1",
		"twigMaterial":"BranchType6"}
		
	);
	
	newTree1.seed=seed[crossSel[0]];
	newTree1.initalBranchLength=initalBranchLength[crossSel[1]];
	newTree1.lengthFalloffFactor=lengthFalloffFactor[crossSel[2]];
	newTree1.lengthFalloffPower=lengthFalloffPower[crossSel[3]];
	newTree1.clumpMax=clumpMax[crossSel[4]];
	newTree1.clumpMin=clumpMin[crossSel[5]];
	newTree1.branchFactor=branchFactor[crossSel[6]];
	newTree1.dropAmount=dropAmount[crossSel[7]];
	newTree1.growAmount=growAmount[crossSel[8]];
	newTree1.sweepAmount=sweepAmount[crossSel[9]];
	newTree1.climbRate=climbRate[crossSel[10]];
	newTree1.trunkKink=trunkKink[crossSel[11]];
	newTree1.taperRate=taperRate[crossSel[12]];
	newTree1.radiusFalloffRate=radiusFalloffRate[crossSel[13]];
	newTree1.twistRate=twistRate[crossSel[14]];
	newTree1.trunkLength=trunkLength[crossSel[15]];
	
	var newTree2 = new Tree(
		{"seed":seed[Math.abs(crossSel[0]-1)],
		"segments":6,
		"levels":5,
		"vMultiplier":1.16,
		"twigScale":0.2,
		"initalBranchLength":initalBranchLength[Math.abs(crossSel[1]-1)],
		"lengthFalloffFactor":lengthFalloffFactor[Math.abs(crossSel[2]-1)],
		"lengthFalloffPower":lengthFalloffPower[Math.abs(crossSel[3]-1)],
		"clumpMax":clumpMax[Math.abs(crossSel[4]-1)],
		"clumpMin":clumpMin[Math.abs(crossSel[5]-1)],
		"branchFactor":branchFactor[Math.abs(crossSel[6]-1)],
		"dropAmount":dropAmount[Math.abs(crossSel[7]-1)],
		"growAmount":growAmount[Math.abs(crossSel[8]-1)],
		"sweepAmount":sweepAmount[Math.abs(crossSel[9]-1)],
		"maxRadius":0.111,
		"climbRate":climbRate[Math.abs(crossSel[10]-1)],
		"trunkKink":trunkKink[Math.abs(crossSel[11]-1)],
		"treeSteps":4,
		"taperRate":taperRate[Math.abs(crossSel[12]-1)],
		"radiusFalloffRate":radiusFalloffRate[Math.abs(crossSel[13]-1)],
		"twistRate":twistRate[Math.abs(crossSel[14]-1)],
		"trunkLength":trunkLength[Math.abs(crossSel[15]-1)],
		"trunkMaterial":"TrunkType1",
		"twigMaterial":"BranchType6"}
		
	);
	
	newTree2.seed=seed[Math.abs(crossSel[0]-1)];
	newTree2.initalBranchLength=initalBranchLength[Math.abs(crossSel[1]-1)];
	newTree2.lengthFalloffFactor=lengthFalloffFactor[Math.abs(crossSel[2]-1)];
	newTree2.lengthFalloffPower=lengthFalloffPower[Math.abs(crossSel[3]-1)];
	newTree2.clumpMax=clumpMax[Math.abs(crossSel[4]-1)];
	newTree2.clumpMin=clumpMin[Math.abs(crossSel[5]-1)];
	newTree2.branchFactor=branchFactor[Math.abs(crossSel[6]-1)];
	newTree2.dropAmount=dropAmount[Math.abs(crossSel[7]-1)];
	newTree2.growAmount=growAmount[Math.abs(crossSel[8]-1)];
	newTree2.sweepAmount=sweepAmount[Math.abs(crossSel[9]-1)];
	newTree2.climbRate=climbRate[Math.abs(crossSel[10]-1)];
	newTree2.trunkKink=trunkKink[Math.abs(crossSel[11]-1)];
	newTree2.taperRate=taperRate[Math.abs(crossSel[12]-1)];
	newTree2.radiusFalloffRate=radiusFalloffRate[Math.abs(crossSel[13]-1)];
	newTree2.twistRate=twistRate[Math.abs(crossSel[14]-1)];
	newTree2.trunkLength=trunkLength[Math.abs(crossSel[15]-1)];
	//console.log("verts: "+myTree.verts.length+", polygons: "+myTree.faces.length+", twigPolygons: "+myTree.facesTwig.length);
	var newtrees=[];
	
	newtrees.push(newTree1);
	newtrees.push(newTree2);
	return newtrees;
}

function crossOver2(tree1,tree2){
	var seed=[tree1.seed, tree2.seed];
	//console.log(seed);
	var initalBranchLength=[tree1.initalBranchLength, tree2.initalBranchLength];
	//console.log(initalBranchLength);
	var lengthFalloffFactor=[tree1.lengthFalloffFactor, tree2.lengthFalloffFactor];
	//console.log(lengthFalloffFactor);
	var lengthFalloffPower=[tree1.lengthFalloffPower, tree2.lengthFalloffPower];
	//console.log(lengthFalloffPower);
	var clumpMax=[tree1.clumpMax, tree2.clumpMax];
	//console.log(clumpMax);
	var clumpMin=[tree1.clumpMin, tree2.clumpMin];
    //console.log(clumpMin);
	var branchFactor=[tree1.branchFactor, tree2.branchFactor];
	//console.log(branchFactor);
	var dropAmount=[tree1.dropAmount, tree2.dropAmount];
	//console.log(dropAmount);
	var growAmount=[tree1.growAmount, tree2.growAmount];
	//console.log(growAmount);
	var sweepAmount=[tree1.sweepAmount, tree2.sweepAmount];
	//console.log(sweepAmount);
	var climbRate=[tree1.climbRate, tree2.climbRate];
	//console.log(climbRate);
	var trunkKink=[tree1.trunkKink, tree2.trunkKink];
	//console.log(trunkKink);
	var taperRate=[tree1.taperRate, tree2.taperRate];
	//console.log(taperRate);
	var radiusFalloffRate=[tree1.radiusFalloffRate, tree2.radiusFalloffRate];
	//console.log(radiusFalloffRate);
	var twistRate=[tree1.twistRate, tree2.twistRate];
	//console.log(twistRate);
	var trunkLength=[tree1.trunkLength, tree2.trunkLength];
	//console.log(trunkLength);
	
	var newTree = new Tree(
		{"seed":Math.round(seed[0]+seed[1])/2,
		"segments":6,
		"levels":5,
		"vMultiplier":1.16,
		"twigScale":0.2,
		"initalBranchLength":(initalBranchLength[0]+initalBranchLength[1])/2,
		"lengthFalloffFactor":(lengthFalloffFactor[0]+lengthFalloffFactor[1])/2,
		"lengthFalloffPower":(lengthFalloffPower[0]+lengthFalloffPower[1])/2,
		"clumpMax":(clumpMax[0]+clumpMax[1])/2,
		"clumpMin":(clumpMin[0]+clumpMin[1])/2,
		"branchFactor":(branchFactor[0]+branchFactor[1])/2,
		"dropAmount":(dropAmount[0]+dropAmount[1])/2,
		"growAmount":(growAmount[0]+growAmount[1])/2,
		"sweepAmount":(sweepAmount[0]+sweepAmount[1])/2,
		"maxRadius":0.111,
		"climbRate":(climbRate[0]+climbRate[1])/2,
		"trunkKink":(trunkKink[0]+trunkKink[1])/2,
		"treeSteps":4,
		"taperRate":(taperRate[0]+taperRate[1])/2,
		"radiusFalloffRate":(radiusFalloffRate[0]+radiusFalloffRate[1])/2,
		"twistRate":(twistRate[0]+twistRate[1])/2,
		"trunkLength":(trunkLength[0]+trunkLength[1])/2,
		"trunkMaterial":"TrunkType1",
		"twigMaterial":"BranchType6"}
		
	);
	
	newTree.seed=(seed[0]+seed[1])/2;
	newTree.initalBranchLength=(initalBranchLength[0]+initalBranchLength[1])/2;
	newTree.lengthFalloffFactor=(lengthFalloffFactor[0]+lengthFalloffFactor[1])/2;
	newTree.lengthFalloffPower=(lengthFalloffPower[0]+lengthFalloffPower[1])/2;
	newTree.clumpMax=(clumpMax[0]+clumpMax[1])/2;
	newTree.clumpMin=(clumpMin[0]+clumpMin[1])/2;
	newTree.branchFactor=(branchFactor[0]+branchFactor[1])/2;
	newTree.dropAmount=(dropAmount[0]+dropAmount[1])/2;
	newTree.growAmount=(growAmount[0]+growAmount[1])/2;
	newTree.sweepAmount=(sweepAmount[0]+sweepAmount[1])/2;
	newTree.climbRate=(climbRate[0]+climbRate[1])/2;
	newTree.trunkKink=(trunkKink[0]+trunkKink[1])/2;
	newTree.taperRate=(taperRate[0]+taperRate[1])/2;
	newTree.radiusFalloffRate=(radiusFalloffRate[0]+radiusFalloffRate[1])/2;
	newTree.twistRate=(twistRate[0]+twistRate[1])/2;
	newTree.trunkLength=(trunkLength[0]+trunkLength[1])/2;
	
	var newtrees=[];
	newtrees.push(newTree);
	return newtrees;
}


function mutateTree(myTree, mutationAmount){
	var rFac=Math.random();

	var seed=Math.round(Math.random()*4000);
	if(rFac>mutationAmount){
		seed=myTree.seed;
	}
	
	rFac=Math.random();
	var initalBranchLength=0.5+Math.random()*0.3;
	if(rFac>mutationAmount){
		initalBranchLength=myTree.initalBranchLength;
	}
		
	rFac=Math.random();
	var lengthFalloffFactor=0.5+Math.random()*0.3;
	if(rFac>mutationAmount){
		lengthFalloffFactor=myTree.lengthFalloffFactor;
	}
		
	rFac=Math.random();
	var lengthFalloffPower=0.3+Math.random()*0.4;
	if(rFac>mutationAmount){
		lengthFalloffPower=myTree.lengthFalloffPower;
	}
		
	rFac=Math.random();
	var clumpMax=0.4+Math.random()*0.1;
	if(rFac>mutationAmount){
		clumpMax=myTree.clumpMax;
	}
		
	rFac=Math.random();
	var clumpMin=clumpMax-Math.random()*0.4;
	if(rFac>mutationAmount){
		clumpMin=myTree.clumpMin;
	}
		
	rFac=Math.random();
	var branchFactor=2.0+Math.random()*2.0;
	if(rFac>mutationAmount){
		branchFactor=myTree.branchFactor;
	}
		
	rFac=Math.random();
	var dropAmount=-0.3+Math.random()*0.6;
	if(rFac>mutationAmount){
		dropAmount=myTree.dropAmount;
	}
		
	rFac=Math.random();
	var growAmount=-0.5+Math.random()*1.5;
	if(rFac>mutationAmount){
		growAmount=myTree.growAmount;
	}
		
	rFac=Math.random();
	var sweepAmount=-0.05+Math.random()*0.1;
	if(rFac>mutationAmount){
		sweepAmount=myTree.sweepAmount;
	}
		
	rFac=Math.random();
	var climbRate=0.05+Math.random()*0.95;
	if(rFac>mutationAmount){
		climbRate=myTree.climbRate;
	}
		
	rFac=Math.random();
	var trunkKink=Math.random()*0.3;
	if(rFac>mutationAmount){
		trunkKink=myTree.trunkKink;
	}
		
	rFac=Math.random();
	var taperRate=0.7+Math.random()*0.3;
	if(rFac>mutationAmount){
		taperRate=myTree.taperRate;
	}
		
	rFac=Math.random();
	var radiusFalloffRate=0.74+Math.random()*0.05;
	if(rFac>mutationAmount){
		radiusFalloffRate=myTree.radiusFalloffRate;
	}
		
	rFac=Math.random();
	var twistRate=Math.random()*10.0;
	if(rFac>mutationAmount){
		twistRate=myTree.twistRate;
	}
		
	rFac=Math.random();
	var trunkLength=1.5+Math.random()*1.6;
	if(rFac>mutationAmount){
		trunkLength=myTree.trunkLength;
	}
	
	var mutatedTree = new Tree(
		{"seed":seed,
		"segments":6,
		"levels":5,
		"vMultiplier":1.16,
		"twigScale":0.22,
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
	
	mutatedTree.seed=seed;
	mutatedTree.initalBranchLength=initalBranchLength;
	mutatedTree.lengthFalloffFactor=lengthFalloffFactor;
	mutatedTree.lengthFalloffPower=lengthFalloffPower;
	mutatedTree.clumpMax=clumpMax;
	mutatedTree.clumpMin=clumpMin;
	mutatedTree.branchFactor=branchFactor;
	mutatedTree.dropAmount=dropAmount;
	mutatedTree.growAmount=growAmount;
	mutatedTree.sweepAmount=sweepAmount;
	mutatedTree.climbRate=climbRate;
	mutatedTree.trunkKink=trunkKink;
	mutatedTree.taperRate=taperRate;
	mutatedTree.radiusFalloffRate=radiusFalloffRate;
	mutatedTree.twistRate=twistRate;
	mutatedTree.trunkLength=trunkLength;
	//console.log("verts: "+myTree.verts.length+", polygons: "+myTree.faces.length+", twigPolygons: "+myTree.facesTwig.length);
	
	return mutatedTree;
}

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
		"twigScale":0.22,
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