//Donald Elmore
//triangle rendering

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define numRows 256
#define numCols 256	
#define PI 3.14159265
double deg2rad = PI/180.0;

typedef struct plyVector {
	double x, y, z;
} vertex;

typedef struct plyTriangle  {
	vertex a, b, c;	//each vertex of the triangle
} triangle;

double dotProduct(vertex a, vertex b) {
	return (a.x * b.x  +  a.y * b.y  +  a.z * b.z);
}

vertex crossProduct(vertex a, vertex b) {
	vertex answer;
	answer.x = a.y * b.z - a.z * b.y;
	answer.y = a.z * b.x - a.x * b.z;
	answer.z = a.x * b.y - a.y * b.x;

	return answer;
}

int main(int argc, char *argv[]) {
	if (argc != 5) {
        printf("Error! Formatting is: ./executable filename.ply xDeg yDeg zDeg\n");
        exit(0);
    }

    int a, b, c, dummy;
	char dum1[100], num[20];
	int numVertices, numFaces, i;
    FILE *in = fopen(argv[1], "rb");

/* Determine number of vertices and faces */
    //eat two lines then read numVertices
    fgets(dum1, 100, in);
    fgets(dum1, 100, in);
    fscanf(in, "%s %s %s", dum1, dum1, num);
    numVertices = atoi(num);

    //eat newline char and 3 more lines then read numFaces
    fgets(dum1, 100, in);
    fgets(dum1, 100, in);
    fgets(dum1, 100, in);
    fgets(dum1, 100, in);
    fscanf(in, "%s %s %s", dum1, dum1, num);
    numFaces = atoi(num);
    fgets(dum1, 100, in);	//eat newline char
    fgets(dum1, 100, in);
    fgets(dum1, 100, in); //eat two last lines before vertices
    //printf("numVertices = %d, numFaces = %d\n", numVertices, numFaces);


/* Read PLY in vertices and faces */
    vertex *vertices = (vertex *)calloc(numVertices, sizeof(vertex));
    triangle *triangles = (triangle *)calloc(numFaces, sizeof(triangle));

    //scan in vertices
    for (i = 0; i < numVertices; i++) {
    	fscanf(in, "%lf %lf %lf", &vertices[i].x, &vertices[i].y, &vertices[i].z);
    	//printf("Vertex %d: %f %f %f\n", i, vertices[i].x, vertices[i].y, vertices[i].z);
    }
    //scan in triangles
    for (i = 0; i < numFaces; i++) {
    	fscanf(in, "%d %d %d %d", &dummy, &a, &b, &c);
    	triangles[i].a = vertices[a];
    	triangles[i].b = vertices[b];
    	triangles[i].c = vertices[c];
    	//printf("triangle %d\n", i);
    	//printf("\tvertex a=%f,%f,%f\n",triangles[i].a.x,triangles[i].a.y,triangles[i].a.z);
    	//printf("\tvertex b=%f,%f,%f\n",triangles[i].b.x,triangles[i].b.y,triangles[i].b.z);
    	//printf("\tvertex c=%f,%f,%f\n",triangles[i].c.x,triangles[i].c.y,triangles[i].c.z);
    }
    fclose(in);
    
/* Calculate the bounding box on the vertices */
	vertex max, min, center;

	//find min vertex, initialize to arb large value so gets reset
	min.x = 99999; min.y = 99999; min.z = 99999;
	for (i = 0; i < numVertices; i++) {
		if (min.x > vertices[i].x)
			min.x =  vertices[i].x;
		if (min.y > vertices[i].y)
			min.y =  vertices[i].y;
		if (min.z > vertices[i].z)
			min.z =  vertices[i].z;
	}
	printf("min = %f,%f,%f\n", min.x, min.y, min.z);

	//find max vertex, initialize to arb small value so gets reset
	max.x = -99999; max.y = -99999; max.z = -99999;
	for (i = 0; i < numVertices; i++) {
		if (max.x < vertices[i].x)
			max.x =  vertices[i].x;
		if (max.y < vertices[i].y)
			max.y =  vertices[i].y;
		if (max.z < vertices[i].z)
			max.z =  vertices[i].z;
	}
	printf("max = %f,%f,%f\n", max.x, max.y, max.z);

	//find center 
	center.x = (max.x + min.x) / 2.0;
	center.y = (max.y + min.y) / 2.0;
	center.z = (max.z + min.z) / 2.0;
	printf("center = %f,%f,%f\n", center.x, center.y, center.z);

	//find maximum extent of bounding box (E)
	double E = fabsf(max.x - min.x);
	if (E < fabsf(max.y - min.y))
		E = fabsf(max.y - min.y);
	if (E < fabsf(max.z - min.z))
		E = fabsf(max.z - min.z);
	printf("E = %f\n", E);


/* Calculate the camera position and orientation using two vectors (camera & up)  */
	vertex camera, up;
	camera.x = 1.0; camera.y = 0.0; camera.z = 0.0;		//default camera position
	up.x = 0.0; up.y = 0.0; up.z = 1.0;					//default up position
	double xDeg = (double)atoi(argv[2]), yDeg = (double)atoi(argv[3]), zDeg = (double)atoi(argv[4]);
	printf("\nInput angles: %f %f %f\n", xDeg, yDeg, zDeg);


	//setup rotation matrices
	double xRot[9], yRot[9], zRot[9];

	xRot[0] = 1.0; 				xRot[1] = 0.0;						xRot[2] = 0.0;
	xRot[3] = 0.0; 				xRot[4] = cosf(xDeg*deg2rad);		xRot[5] = -1.0*sinf(xDeg*deg2rad);
	xRot[6] = 0.0; 				xRot[7] = sinf(xDeg*deg2rad); 		xRot[8] = cosf(xDeg*deg2rad);

	yRot[0] = cosf(yDeg*deg2rad); 		yRot[1] = 0.0;				yRot[2] = sinf(yDeg*deg2rad);
	yRot[3] = 0.0; 						yRot[4] = 1.0;				yRot[5] = 0.0;
	yRot[6] = -1.0*sinf(yDeg*deg2rad); 	yRot[7] = 0.0; 				yRot[8] = cosf(yDeg*deg2rad);

	zRot[0] = cosf(zDeg*deg2rad); 		zRot[1] = -1.0*sinf(zDeg*deg2rad);	zRot[2] = 0.0;
	zRot[3] = sinf(zDeg*deg2rad); 		zRot[4] = cosf(zDeg*deg2rad);		zRot[5] = 0.0;
	zRot[6] = 0.0; 						zRot[7] = 0.0; 						zRot[8] = 1.0;

	//rotate camera and up vectors by rotations entered at command line
	printf("\nRotations:\n");
	vertex ta;
	ta.x = camera.x; ta.y = camera.y; ta.z = camera.z;
	camera.x = ta.x * xRot[0] + ta.y * xRot[1] + ta.z * xRot[2];
	camera.y = ta.x * xRot[3] + ta.y * xRot[4] + ta.z * xRot[5];
	camera.z = ta.x * xRot[6] + ta.y * xRot[7] + ta.z * xRot[8];
	printf("camera1st: %f,%f,%f\n", camera.x,camera.y,camera.z);

	ta.x = camera.x; ta.y = camera.y; ta.z = camera.z;
	camera.x = ta.x * yRot[0] + ta.y * yRot[1] + ta.z * yRot[2];
	camera.y = ta.x * yRot[3] + ta.y * yRot[4] + ta.z * yRot[5];
	camera.z = ta.x * yRot[6] + ta.y * yRot[7] + ta.z * yRot[8];
	printf("camera2st: %f,%f,%f\n", camera.x,camera.y,camera.z);

	ta.x = camera.x; ta.y = camera.y; ta.z = camera.z;
	camera.x = ta.x * zRot[0] + ta.y * zRot[1] + ta.z * zRot[2];
	camera.y = ta.x * zRot[3] + ta.y * zRot[4] + ta.z * zRot[5];
	camera.z = ta.x * zRot[6] + ta.y * zRot[7] + ta.z * zRot[8];
	printf("camera3st: %f,%f,%f\n", camera.x,camera.y,camera.z);

	//up
	ta.x = up.x; ta.y = up.y; ta.z = up.z;
	up.x = ta.x * xRot[0] + ta.y * xRot[1] + ta.z * xRot[2];
	up.y = ta.x * xRot[3] + ta.y * xRot[4] + ta.z * xRot[5];
	up.z = ta.x * xRot[6] + ta.y * xRot[7] + ta.z * xRot[8];
	printf("up1st: %f,%f,%f\n", up.x,up.y,up.z);

	ta.x = up.x; ta.y = up.y; ta.z = up.z;
	up.x = ta.x * yRot[0] + ta.y * yRot[1] + ta.z * yRot[2];
	up.y = ta.x * yRot[3] + ta.y * yRot[4] + ta.z * yRot[5];
	up.z = ta.x * yRot[6] + ta.y * yRot[7] + ta.z * yRot[8];
	printf("up2st: %f,%f,%f\n", up.x,up.y,up.z);

	ta.x = up.x; ta.y = up.y; ta.z = up.z;
	up.x = ta.x * zRot[0] + ta.y * zRot[1] + ta.z * zRot[2];
	up.y = ta.x * zRot[3] + ta.y * zRot[4] + ta.z * zRot[5];
	up.z = ta.x * zRot[6] + ta.y * zRot[7] + ta.z * zRot[8];
	printf("up3st: %f,%f,%f\n\n", up.x,up.y,up.z);

	//scale camera vector
	camera.x = 1.5 * E * camera.x + center.x;
	camera.y = 1.5 * E * camera.y + center.y;
	camera.z = 1.5 * E * camera.z + center.z;
	printf("camera: %f,%f,%f\n", camera.x,camera.y,camera.z);


/* Determine the 3D coordinates bounding the image */
	double A;
	vertex left, right, bottom, top, topleft, temp, temp2, temp3;
	memset(&left, 0, sizeof(left));
	memset(&right, 0, sizeof(left));
	memset(&bottom, 0, sizeof(left));
	memset(&top, 0, sizeof(left));
	memset(&topleft, 0, sizeof(left));
	memset(&temp, 0, sizeof(left));

	//left = up x (center - camera)
	temp.x = center.x - camera.x;
	temp.y = center.y - camera.y;
	temp.z = center.z - camera.z;
	printf("temp = %f,%f,%f\n", temp.x, temp.y, temp.z);
	left = crossProduct(up, temp);
	printf("left1 = %f,%f,%f\n", left.x, left.y, left.z);

	//A = ||left||
	A = sqrtf(left.x * left.x + left.y * left.y + left.z * left.z);
	printf("A = %f\n", A);

	//left = E/2A * left + center
	left.x = ((E / (2.0 * A)) * left.x) + center.x;
	left.y = ((E / (2.0 * A)) * left.y) + center.y;
	left.z = ((E / (2.0 * A)) * left.z) + center.z;
	printf("left2 = %f,%f,%f\n", left.x, left.y, left.z);
 
 	//right = (center - camera) x up... reusing temp from above as (center - camera)
	right = crossProduct(temp, up);
	printf("right1 = %f,%f,%f\n", right.x, right.y, right.z);

	//right = E/2A * right + center
	right.x = ((E / (2.0 * A)) * right.x) + center.x;
	right.y = ((E / (2.0 * A)) * right.y) + center.y;
	right.z = ((E / (2.0 * A)) * right.z) + center.z;
	printf("right2 = %f,%f,%f\n", right.x, right.y, right.z);

	//top = E/2 * up + center
	top.x = ((E / 2.0) * up.x) + center.x;
	top.y = ((E / 2.0) * up.y) + center.y;
	top.z = ((E / 2.0) * up.z) + center.z;
	printf("top = %f,%f,%f\n", top.x, top.y, top.z);

	//bottom = -E/2 * up + center
	bottom.x = ((-1.0 * E / 2.0) * up.x) + center.x;
	bottom.y = ((-1.0 * E / 2.0) * up.y) + center.y;
	bottom.z = ((-1.0 * E / 2.0) * up.z) + center.z;
	printf("bottom = %f,%f,%f\n", bottom.x, bottom.y, bottom.z);

	//topleft = E/2 * up + left
	topleft.x = ((E / 2.0) * up.x) + left.x;
	topleft.y = ((E / 2.0) * up.y) + left.y;
	topleft.z = ((E / 2.0) * up.z) + left.z;
	printf("topleft = %f,%f,%f\n", topleft.x, topleft.y, topleft.z);


/* For each pixel r,c in the image calculate needed values */
	int R, C, T;											//image size
	unsigned char *image = (unsigned char *)calloc(65536, sizeof(unsigned char));			//default image color is black
	for (int init = 0; init < 65536; init++)
		image[init] = (unsigned char)0;

	vertex imagePx, plEq, intersect;
	double zBuffer, D, d, n, distance;	//zBuffer depth is far
	memset(&imagePx, 0, sizeof(imagePx));
	memset(&plEq, 0, sizeof(plEq));
	memset(&intersect, 0, sizeof(intersect));

	for (R = 0; R < numRows; R++) {
		for (C = 0; C < numCols; C++) {
			zBuffer = 999999.0;

			//imagePx = topleft + ((c/(numCols-1))*(right-left)) + ((r/(numRows-1))*(bottom-top))
			imagePx.x = topleft.x + ((((double)C)/(((double)numCols)-1))*(right.x-left.x)) + (( ((double)R) /(((double)numRows)-1))*(bottom.x-top.x));
			imagePx.y = topleft.y + ((((double)C)/(((double)numCols)-1))*(right.y-left.y)) + (( ((double)R) /(((double)numRows)-1))*(bottom.y-top.y));
			imagePx.z = topleft.z + ((((double)C)/(((double)numCols)-1))*(right.z-left.z)) + (( ((double)R) /(((double)numRows)-1))*(bottom.z-top.z));
			//printf("image: %f,%f,%f\n", imagePx.x, imagePx.y, imagePx.z);

			//for each triangle
			for (T = 0; T < numFaces; T++) {
				//<A, B, C> = <v1-v0> x <v2-v0>
				temp.x = triangles[T].b.x - triangles[T].a.x;
				temp.y = triangles[T].b.y - triangles[T].a.y;
				temp.z = triangles[T].b.z - triangles[T].a.z;

				temp2.x = triangles[T].c.x - triangles[T].a.x;
				temp2.y = triangles[T].c.y - triangles[T].a.y;
				temp2.z = triangles[T].c.z - triangles[T].a.z;

				plEq = crossProduct(temp, temp2);	//plEq = <A, B, C>
				//printf("plEq = %f,%f,%f\n", plEq.x, plEq.y, plEq.z);

				//D = -<A, B, C> dot <v0>
				temp.x = -1.0 * plEq.x;
				temp.y = -1.0 * plEq.y;
				temp.z = -1.0 * plEq.z;
				D = dotProduct(temp, triangles[T].a);
				//printf("D = %f\n", D);

				//n = -<A, B, C> dot <camera> - D
				n = dotProduct(temp, camera) - D;
				//printf("n = %f\n", n);

				//d = <A, B, C> dot <image - camera>
				temp.x = imagePx.x - camera.x;
				temp.y = imagePx.y - camera.y;
				temp.z = imagePx.z - camera.z;
				d = dotProduct(plEq, temp);
				//printf("d = %f\n", d);

				//check if d is near zero, if so, skip this triangle for this pixel
				//printf("fabs(d) = %f\n", fabs(d));
				if (fabs(d) <= 0.0001) {
					//printf("d is near zero\n");
					continue;
				}
				distance = n/d;
				//printf("distance = %f\n", distance);

				//find 3D coordinates <intersect> of ray and plane
				intersect.x = camera.x + (distance * (imagePx.x - camera.x));
				intersect.y = camera.y + (distance * (imagePx.y - camera.y)); 
				intersect.z = camera.z + (distance * (imagePx.z - camera.z));
				//printf("intersect = %f,%f,%f\n", intersect.x, intersect.y, intersect.z);

				//Determine if intersection point lies within triangle
				double dot1 = 0.0, dot2 = 0.0, dot3 = 0.0;

				//dot1 = <v2-v0> x <v1-v0> dot <intersect-v0> x <v1-v0>
					//<v2-v0>
					temp.x = triangles[T].c.x - triangles[T].a.x;
					temp.y = triangles[T].c.y - triangles[T].a.y;
					temp.z = triangles[T].c.z - triangles[T].a.z;

					//<v1-v0>
					temp2.x = triangles[T].b.x - triangles[T].a.x;
					temp2.y = triangles[T].b.y - triangles[T].a.y;
					temp2.z = triangles[T].b.z - triangles[T].a.z;

					//<intersect-v0>
					temp3.x = intersect.x - triangles[T].a.x;
					temp3.y = intersect.y - triangles[T].a.y;
					temp3.z = intersect.z - triangles[T].a.z;

				dot1 = dotProduct(crossProduct(temp, temp2), crossProduct(temp3, temp2));
				//printf("dot1 = %lf\n", dot1);
				if (dot1 < 0)
					continue;

				//dot2 = <v0-v1> x <v2-v1> dot <intersect-v1> x <v2-v1>
					//<v0-v1>
					temp.x = triangles[T].a.x - triangles[T].b.x;
					temp.y = triangles[T].a.y - triangles[T].b.y;
					temp.z = triangles[T].a.z - triangles[T].b.z;

					//<v2-v1>
					temp2.x = triangles[T].c.x - triangles[T].b.x;
					temp2.y = triangles[T].c.y - triangles[T].b.y;
					temp2.z = triangles[T].c.z - triangles[T].b.z;

					//<intersect-v1>
					temp3.x = intersect.x - triangles[T].b.x;
					temp3.y = intersect.y - triangles[T].b.y;
					temp3.z = intersect.z - triangles[T].b.z;

				dot2 = dotProduct(crossProduct(temp, temp2), crossProduct(temp3, temp2));
				//printf("dot2 = %lf\n", dot2);
				if (dot2 < 0)
					continue;

				//dot3 = <v1-v2> x <v0-v2> dot <intersect-v2> x <v0-v2>
					//<v1-v2>
					temp.x = triangles[T].b.x - triangles[T].c.x;
					temp.y = triangles[T].b.y - triangles[T].c.y;
					temp.z = triangles[T].b.z - triangles[T].c.z;

					//<v0-v2>
					temp2.x = triangles[T].a.x - triangles[T].c.x;
					temp2.y = triangles[T].a.y - triangles[T].c.y;
					temp2.z = triangles[T].a.z - triangles[T].c.z;

					//<intersect-v2>
					temp3.x = intersect.x - triangles[T].c.x;
					temp3.y = intersect.y - triangles[T].c.y;
					temp3.z = intersect.z - triangles[T].c.z;

				dot3 = dotProduct(crossProduct(temp, temp2), crossProduct(temp3, temp2));
				//printf("dot3 = %lf\n", dot3);
				if (dot3 < 0)
					continue;

				//if ((dot1 < 0) || (dot2 < 0) || (dot3 < 0)) {	
					//printf("skip triangle dot < 0\n");
				//	continue;
				//}

				//Check if distance to triangle is greater than current z-buffer
				//printf("distance = %f", distance);
				if (distance > zBuffer) {
					//printf("distance > zBuffer\n");
					continue;
				}
				if (distance < zBuffer) {
					//printf("set pixel value to :%u\n", (unsigned char)(155 + T % 100));
					zBuffer = distance;
					image[numCols * R + C] = (unsigned char)(155 + T % 100);
					//printf("T = %d, image[%d] = %u\n", T, numCols*R+C, image[numCols*R+C]);
				}
			}
			//printf("image1[%d] = %u\n", numCols*R+C, image[numCols*R+C]);
		}
	}

	printf("Print Image: \n");
	for (int p = 0; p < 65536; p++)
		printf("%u", image[p]);
	printf("\nDONE PRINTING\n");

	/* Write PPM Image */
	FILE *outt = fopen("output.ppm", "wb");
	//char h[] = "P5 256 256 255\n";
	fprintf(outt, "P5 256 256 255\n");

	//fwrite(&image, 65536 * sizeof(unsigned char), 1, outt);

	for (int ima; ima < 65536; ima++) {
		fwrite(&image[ima], sizeof(unsigned char), 1, outt);
		printf("%u", image[ima]);
	}

	printf("\n");
    fclose(outt);
    free(image);
	return 0;
}
