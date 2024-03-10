#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <unistd.h>
#include <SDL2/SDL.h>
#include <signal.h>

#define ORANGE 0xFFFFA500 // (aa-bb-gg-rr)
#define RED 0xFFFF0000
#define WHITE 0xFFFFFFFF
#define YELLOW 0xFFFFFF00
#define BLUE 0xFF0000FF
#define GREEN 0xFF00DD00
#define BLACK 0xFF040404
#define BACKGROUND_COLOR 0xFF7FBFFF
#define BACKGROUND_COLOR_DARK 0xFF101820
#define FLOOR_GRID_COLOR 0xFF206020

#define MAX_TRIS 10000

#define WINDOW_WIDTH 600
#define WINDOW_HEIGHT 400
#define PI 3.14159265358979323846

bool drawShadows = false;


// Data structures
typedef struct {
    double x, y, z;
} Vertex;

typedef struct {
    Vertex **tri_address;
    bool *is_quad_tri;
    uint32_t **object_address;
    Vertex **v;
} CollisionTris;

CollisionTris collisionTris;
int collision_n = 0;

typedef struct {
    double m[4][4];
} Matrix;

Matrix matrix_multiply(Matrix a, Matrix b) {
    Matrix result = {{{0}}};
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                result.m[i][j] += a.m[i][k] * b.m[k][j];
            }
        }
    }

    return result;
}


#define fov ((double)1.1)
#define aspect ((double)WINDOW_WIDTH/(double)WINDOW_HEIGHT)
#define near ((double)0.6)
#define far ((double)100.0)
#define tan_fov ((double)0.61310521328) // tanf(1.1 / 2)

// perspective projection matrix
const Matrix projection = {
    .m = {
        {1.0f / (tan_fov * aspect), 0, 0, 0},
        {0, 1.0f / tan_fov, 0, 0},
        {0, 0, -(far + near) / (far - near), -(2 * far * near) / (far - near)},
        {0, 0, -1, 0}
    }
};

// inverse projection matrix
const Matrix deprojection = {
    .m = {
        {tan_fov * aspect, 0, 0, 0},
        {0, tan_fov, 0, 0},
        {0, 0, 0, -1},
        {0, 0, (-1.0f / (2 * far * near)) * (far - near), (far + near) / (2 * far * near)}
    }
};



Vertex vertex_project(Vertex v, Matrix m) {
    Vertex result;
    result.x = (v.x * m.m[0][0]) + (v.y * m.m[1][0]) + (v.z * m.m[2][0]) + m.m[3][0];
    result.y = (v.x * m.m[0][1]) + (v.y * m.m[1][1]) + (v.z * m.m[2][1]) + m.m[3][1];
    result.z = (v.x * m.m[0][2]) + (v.y * m.m[1][2]) + (v.z * m.m[2][2]) + m.m[3][2];

    double w = (v.x * m.m[0][3]) + (v.y * m.m[1][3]) + (v.z * m.m[2][3]) + m.m[3][3];

    if(w != 0) {
        result.x /= w;
        result.y /= w;
        result.z /= w;
    }
    
    return result;
}

double vz_to_pz(double vz) {
    double pz = ((1.0/vz) - projection.m[2][2]) / -projection.m[2][3];
    return pz;
}

double vz_to_pz_with_scalar(double vz, double scalar) {
    double pz = ((1.0/(vz*scalar)) - projection.m[2][2]) / -projection.m[2][3];
    return pz;
}

double pz_to_vz(double pz) {
    double vz = 1.0/((-projection.m[2][3] * pz) + projection.m[2][2]);
    return vz;
}

Vertex vertex_rotate(Vertex v, Matrix m) {
    Vertex result;
    result.x = (v.x * m.m[0][0]) + (v.y * m.m[1][0]) + (v.z * m.m[2][0]);
    result.y = (v.x * m.m[0][1]) + (v.y * m.m[1][1]) + (v.z * m.m[2][1]);
    result.z = (v.x * m.m[0][2]) + (v.y * m.m[1][2]) + (v.z * m.m[2][2]);
    return result;
}

Matrix rotation_axis(Vertex axis, double angle) {
    double c = cos(angle);
    double s = sin(angle);
    double t = 1.0f - c;
    double x = axis.x, y = axis.y, z = axis.z;

    // Normalize the axis vector
    double magnitude = sqrt(x*x + y*y + z*z);
    x /= magnitude;
    y /= magnitude;
    z /= magnitude;

    Matrix rot = {
        .m = {
            {t*x*x + c,     t*x*y - s*z,   t*x*z + s*y,   0},
            {t*x*y + s*z,   t*y*y + c,     t*y*z - s*x,   0},
            {t*x*z - s*y,   t*y*z + s*x,   t*z*z + c,     0},
            {0,             0,             0,             1}
        }
    };
    return rot;
}



Vertex vertex_add(Vertex v1, Vertex v2) {
    Vertex result;
    result.x = v1.x + v2.x;
    result.y = v1.y + v2.y;
    result.z = v1.z + v2.z;
    return result;
}

Vertex vertex_subtract(Vertex v1, Vertex v2) {
    Vertex result;
    result.x = v1.x - v2.x;
    result.y = v1.y - v2.y;
    result.z = v1.z - v2.z;
    return result;
}

Vertex vertex_scale(Vertex v, double scaleFactor) {
    Vertex result;
    result.x = v.x * scaleFactor;
    result.y = v.y * scaleFactor;
    result.z = v.z * scaleFactor;
    return result;
}




// rotation matrices
Matrix rotation_x(double angle) {
    double c = cos(angle);
    double s = sin(angle);
    Matrix rotX = {
        .m = {
            {1, 0, 0, 0},
            {0, c, -s, 0},
            {0, s, c, 0},
            {0, 0, 0, 1}
        }
    };
    return rotX;
}

Matrix rotation_y(double angle) {
    double c = cos(angle);
    double s = sin(angle);
    Matrix rotY = {
        .m = {
            {c, 0, s, 0},
            {0, 1, 0, 0},
            {-s, 0, c, 0},
            {0, 0, 0, 1}
        }
    };
    return rotY;
}

Matrix rotation_z(double angle) {
    double c = cos(angle);
    double s = sin(angle);
    Matrix rotZ = {
        .m = {
            {c, -s, 0, 0},
            {s, c, 0, 0},
            {0, 0, 1, 0},
            {0, 0, 0, 1}
        }
    };
    return rotZ;
}















double radians_to_degrees(double radians) {
    double degrees = radians * (180 / PI);
    return degrees;
}

double degrees_to_radians(double degrees) {
    double radians = degrees * (PI / 180);
    return radians;
}

double fix_angle(double radians, double bias) {
    double degrees = radians_to_degrees(radians);
    // range = [0, 360)
    if(degrees < 0) {
        degrees += 360;
    }
    else if(degrees >= 360) {
        degrees -= 360;
    }

    // snapping if angle is close to cardinal
    if(degrees < 0 + bias || degrees > 360 - bias) {
        degrees = 0;
    }
    else if(degrees < 90 + bias && degrees > 90 - bias) {
        degrees = 90;
    }
    else if(degrees < 180 + bias && degrees > 180 - bias) {
        degrees = 180;
    }
    else if(degrees < 270 + bias && degrees > 270 - bias) {
        degrees = 270;
    }

    radians = degrees_to_radians(degrees);
    return radians;
}

bool is_fixed_90(double radians, double bias) {
    bool fixed = false;
    double degrees = radians_to_degrees(radians);
    // range = [0, 360)
    if(degrees < 0) {
        while(degrees < 0) {
            degrees += 360;
        }
    }
    else if(degrees >= 360) {
        while(degrees >= 360) {
            degrees -= 360;
        }
    }

    // snapping if angle is close to cardinal
    if(degrees < 0 + bias || degrees > 360 - bias) {
        degrees = 0;
        fixed = true;
    }
    else if(degrees < 90 + bias && degrees > 90 - bias) {
        degrees = 90;
        fixed = true;
    }
    else if(degrees < 180 + bias && degrees > 180 - bias) {
        degrees = 180;
        fixed = true;
    }
    else if(degrees < 270 + bias && degrees > 270 - bias) {
        degrees = 270;
        fixed = true;
    }

    return fixed;
}




bool is_fixed_120(double radians, double bias) {
    bool fixed = false;
    double degrees = radians_to_degrees(radians);
    // range = [0, 360)
    if(degrees < 0) {
        while(degrees < 0) {
            degrees += 360;
        }
    }
    else if(degrees >= 360) {
        while(degrees >= 360) {
            degrees -= 360;
        }
    }

    // snapping if angle is close to cardinal
    if(degrees < 0 + bias || degrees > 360 - bias) {
        degrees = 0;
        fixed = true;
    }
    else if(degrees < 120 + bias && degrees > 120 - bias) {
        degrees = 90;
        fixed = true;
    }
    else if(degrees < 240 + bias && degrees > 240 - bias) {
        degrees = 180;
        fixed = true;
    }

    return fixed;
}

























// clipping functions
int clip_vertex_near(Vertex* v1, Vertex* v2) {
        // Ensure that at least one vertex is in front of the near plane
    if (v1->z > -near && v2->z > -near) {
        return 0; // Both vertices are behind the near plane, don’t draw the line
    }
    else if (v1->z > -near) {
        // v1 is behind the near plane, v2 is in front
        double t = (-near - v2->z) / (v1->z - v2->z);
        v1->x = v2->x + t * (v1->x - v2->x);
        v1->y = v2->y + t * (v1->y - v2->y);
        v1->z = -near;
        return 1;
    } 
    else if (v2->z > -near) {
        // v2 is behind the near plane, v1 is in front
        double t = (-near - v1->z) / (v2->z - v1->z);
        v2->x = v1->x + t * (v2->x - v1->x);
        v2->y = v1->y + t * (v2->y - v1->y);
        v2->z = -near;
        return 2;
    }

    return 3; // Indicate that the line should be drawn
}

int clip_vertex_far(Vertex* v2, Vertex* v1) {
    // Ensure that at least one vertex is in front of the far plane
    if (v1->z < -far && v2->z < -far) {
        return 0; // Both vertices are beyond the far plane, don’t draw the line
    } 
    else if (v1->z < -far) {
        // v1 is beyond the far plane, v2 is in front
        double t = (-far - v2->z) / (v1->z - v2->z);
        v1->x = v2->x + t * (v1->x - v2->x);
        v1->y = v2->y + t * (v1->y - v2->y);
        v1->z = -far;
        return 2;
    } 
    else if (v2->z < -far) {
        // v2 is beyond the far plane, v1 is in front
        double t = (-far - v1->z) / (v2->z - v1->z);
        v2->x = v1->x + t * (v2->x - v1->x);
        v2->y = v1->y + t * (v2->y - v1->y);
        v2->z = -far;
        return 1;
    }

    return 3; // Indicate that the line should be drawn
}









// view transformations

// vertex
double y_offset = 0.5;
double x_offset = 0.5;
void view_transform_vertex(Vertex v, SDL_Point *point) {
    (*point).x = (Sint16)((v.x + 1.0) * (double)WINDOW_WIDTH / 2.0);
    (*point).y = (Sint16)((-v.y + 1.0) * (double)WINDOW_HEIGHT / 2.0);
}

// vertex (inverse)
void inv_view_transform_vertex(SDL_Point point, double *px, double *py) {
    *px = (((double)point.x + x_offset) * 2.0 / (double)WINDOW_WIDTH) - 1.0;
    *py = -((((double)point.y + y_offset) * 2.0 / (double)WINDOW_HEIGHT) - 1.0);
}

void inv_view_transform_vertex_f(int x, int y, double *px, double *py) {
    *px = (((double)x * 2.0) / (double)WINDOW_WIDTH) - 1.0;
    *py = -((((double)y * 2.0) / (double)WINDOW_HEIGHT) - 1.0);
}

// line
void view_transform_line(Vertex v[2], SDL_Point points[2]) {
    for(int i = 0; i < 2; i++) {
        points[i].x = (Sint16)((v[i].x + 1.0) * (double)WINDOW_WIDTH / 2.0);
        points[i].y = (Sint16)((-v[i].y + 1.0) * (double)WINDOW_HEIGHT / 2.0);
    }
}

// tri
void view_transform_tri(Vertex v[3], SDL_Point points[3]) {
    for(int i = 0; i < 3; i++) {
        points[i].x = (Sint16)((v[i].x + 1.0) * (double)WINDOW_WIDTH / 2.0);
        points[i].y = (Sint16)((-v[i].y + 1.0) * (double)WINDOW_HEIGHT / 2.0);
    }
}



Vertex center_of_polygon(Vertex v[], int n) {
    double x;
    double y;
    double z;
    for(int i = 0; i < n; i++) {
        x += v[i].x;
        y += v[i].y;
        z += v[i].z;
    }
    x /= (double)n;
    y /= (double)n;
    z /= (double)n;
    Vertex center = {x, y, z};
    return center;
}

Vertex center_of_tri(Vertex v[3]) {
    double x = (v[0].x + v[1].x + v[2].x) / 3.0;
    double y = (v[0].y + v[1].y + v[2].y) / 3.0;
    double z = (v[0].z + v[1].z + v[2].z) / 3.0;
    Vertex center = {x, y, z};
    return center;
}



































typedef struct {
    Vertex points[2];
} Line;



// cube stuff
typedef struct {
    Vertex tri[3];
    int color;
    bool isCounterClockwise;
} Face;

typedef struct {
    Vertex vertices[8];
    Vertex edges[12][2];
    Face faces[12];
} Cube;

Cube createCube(Vertex v, double s) {
    
    Cube cube = {
        .vertices = {
            {-0.5 * s + v.x, -0.5 * s + v.y, -0.5 * s + v.z},
            {0.5 * s + v.x, -0.5 * s + v.y, -0.5 * s + v.z},
            {0.5 * s + v.x, 0.5 * s + v.y, -0.5 * s + v.z},
            {-0.5 * s + v.x, 0.5 * s + v.y, -0.5 * s + v.z},
            {-0.5 * s + v.x, -0.5 * s + v.y, 0.5 * s + v.z},
            {0.5 * s + v.x, -0.5 * s + v.y, 0.5 * s + v.z},
            {0.5 * s + v.x, 0.5 * s + v.y, 0.5 * s + v.z},
            {-0.5 * s + v.x, 0.5 * s + v.y, 0.5 * s + v.z}
        },
        .edges = {
            {cube.vertices[0], cube.vertices[1]}, // 0, 1
            {cube.vertices[1], cube.vertices[2]}, // 1, 2
            {cube.vertices[2], cube.vertices[3]}, // 2, 3
            {cube.vertices[3], cube.vertices[0]}, // 3, 0
            {cube.vertices[4], cube.vertices[5]}, // 4, 5
            {cube.vertices[5], cube.vertices[6]}, // 5, 6
            {cube.vertices[6], cube.vertices[7]}, // 6, 7
            {cube.vertices[7], cube.vertices[4]}, // 7, 4
            {cube.vertices[0], cube.vertices[4]}, // 0, 4
            {cube.vertices[1], cube.vertices[5]}, // 1, 5
            {cube.vertices[2], cube.vertices[6]}, // 2, 6
            {cube.vertices[3], cube.vertices[7]} // 3, 7
        },
        .faces[0] = {   
            .tri = {cube.vertices[0], cube.vertices[1], cube.vertices[2]}, // front face
            .color = ORANGE
        },
        .faces[1] = {   
            .tri = {cube.vertices[2], cube.vertices[3], cube.vertices[0]}, // front face
            .color = ORANGE
        },

        .faces[2] = {
            .tri = {cube.vertices[6], cube.vertices[5], cube.vertices[4]}, // back face
            .color = RED
        },
        .faces[3] = {
            .tri = {cube.vertices[4], cube.vertices[7], cube.vertices[6]}, // back face
            .color = RED
        },

        .faces[4] = {
            .tri = {cube.vertices[5], cube.vertices[1], cube.vertices[0]}, // bottom face
            .color = WHITE
        },
        .faces[5] = {
            .tri = {cube.vertices[0], cube.vertices[4], cube.vertices[5]}, // bottom face
            .color = WHITE
        },

        .faces[6] = {
            .tri = {cube.vertices[3], cube.vertices[2], cube.vertices[6]}, // top face
            .color = YELLOW
        },
        .faces[7] = {
            .tri = {cube.vertices[6], cube.vertices[7], cube.vertices[3]}, // top face
            .color = YELLOW
        },

        .faces[8] = {
            .tri = {cube.vertices[0], cube.vertices[3], cube.vertices[7]}, // left face
            .color = BLUE
        },
        .faces[9] = {
            .tri = {cube.vertices[7], cube.vertices[4], cube.vertices[0]}, // left face
            .color = BLUE
        },

        .faces[10] = {
            .tri = {cube.vertices[6], cube.vertices[2], cube.vertices[1]}, // right face
            .color = GREEN
        },
        .faces[11] = {
            .tri = {cube.vertices[1], cube.vertices[5], cube.vertices[6]}, // right face
            .color = GREEN
        },
    };
    return cube;
}


Cube updateCubeEdgesAndFaces(Cube cube) {
    Cube new_cube = {
        .vertices = {
            cube.vertices[0],
            cube.vertices[1],
            cube.vertices[2],
            cube.vertices[3],
            cube.vertices[4],
            cube.vertices[5],
            cube.vertices[6],
            cube.vertices[7]
        },
        .edges = {
            {cube.vertices[0], cube.vertices[1]}, // 0, 1
            {cube.vertices[1], cube.vertices[2]}, // 1, 2
            {cube.vertices[2], cube.vertices[3]}, // 2, 3
            {cube.vertices[3], cube.vertices[0]}, // 3, 0
            {cube.vertices[4], cube.vertices[5]}, // 4, 5
            {cube.vertices[5], cube.vertices[6]}, // 5, 6
            {cube.vertices[6], cube.vertices[7]}, // 6, 7
            {cube.vertices[7], cube.vertices[4]}, // 7, 4
            {cube.vertices[0], cube.vertices[4]}, // 0, 4
            {cube.vertices[1], cube.vertices[5]}, // 1, 5
            {cube.vertices[2], cube.vertices[6]}, // 2, 6
            {cube.vertices[3], cube.vertices[7]} // 3, 7
        },
        .faces[0] = {   
            .tri = {cube.vertices[0], cube.vertices[1], cube.vertices[2]}, // front face
            .color = ORANGE
        },
        .faces[1] = {   
            .tri = {cube.vertices[2], cube.vertices[3], cube.vertices[0]}, // front face
            .color = ORANGE
        },

        .faces[2] = {
            .tri = {cube.vertices[6], cube.vertices[5], cube.vertices[4]}, // back face
            .color = RED
        },
        .faces[3] = {
            .tri = {cube.vertices[4], cube.vertices[7], cube.vertices[6]}, // back face
            .color = RED
        },

        .faces[4] = {
            .tri = {cube.vertices[5], cube.vertices[1], cube.vertices[0]}, // bottom face
            .color = WHITE
        },
        .faces[5] = {
            .tri = {cube.vertices[0], cube.vertices[4], cube.vertices[5]}, // bottom face
            .color = WHITE
        },

        .faces[6] = {
            .tri = {cube.vertices[3], cube.vertices[2], cube.vertices[6]}, // top face
            .color = YELLOW
        },
        .faces[7] = {
            .tri = {cube.vertices[6], cube.vertices[7], cube.vertices[3]}, // top face
            .color = YELLOW
        },

        .faces[8] = {
            .tri = {cube.vertices[0], cube.vertices[3], cube.vertices[7]}, // left face
            .color = BLUE
        },
        .faces[9] = {
            .tri = {cube.vertices[7], cube.vertices[4], cube.vertices[0]}, // left face
            .color = BLUE
        },

        .faces[10] = {
            .tri = {cube.vertices[6], cube.vertices[2], cube.vertices[1]}, // right face
            .color = GREEN
        },
        .faces[11] = {
            .tri = {cube.vertices[1], cube.vertices[5], cube.vertices[6]}, // right face
            .color = GREEN
        },
    };
    return new_cube;
}



















#define RUBIKS_CUBE_DIM 3

typedef struct {
    Cube ***cubies;
    double layer_angle[3][RUBIKS_CUBE_DIM];

} RubiksCube;

RubiksCube createRubiksCube(RubiksCube rubiksCube, double scale) {
    rubiksCube.cubies = malloc((RUBIKS_CUBE_DIM) * sizeof(Cube**));
    for(int i = 0; i < RUBIKS_CUBE_DIM; i++) {
        rubiksCube.cubies[i] = malloc((RUBIKS_CUBE_DIM) * sizeof(Cube*));
        for(int j = 0; j < RUBIKS_CUBE_DIM; j++) {
            rubiksCube.cubies[i][j] = malloc((RUBIKS_CUBE_DIM) * sizeof(Cube));
        }
    }
    double offset = ((double)RUBIKS_CUBE_DIM / 2)*scale - (scale/2); // to maintain center at 0,0,0
    Vertex v;
    for(int i = 0; i < RUBIKS_CUBE_DIM; i++) {
        for(int j = 0; j < RUBIKS_CUBE_DIM; j++) {
            for(int k = 0; k < RUBIKS_CUBE_DIM; k++) {
                v.x = (double)i*scale - offset; 
                v.y = (double)j*scale - offset;
                v.z = (double)k*scale - offset;
                rubiksCube.cubies[i][j][k] = createCube(v, scale);
            }
        }
    }

    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < RUBIKS_CUBE_DIM; j++) {
            rubiksCube.layer_angle[i][j] = 0.0f;
        }
    }



    return rubiksCube;
}

RubiksCube rotateRubiksCube(RubiksCube rubiksCube, double turn_spd, int axis, int layer, int direction) {

    turn_spd *= direction;
    rubiksCube.layer_angle[axis][layer] += turn_spd;
    if(axis == 0) {
        Matrix rotation = rotation_x(turn_spd);
        for(int i = layer; i < (layer + 1); i++) {
            for(int j = 0; j < RUBIKS_CUBE_DIM; j++) {
                for(int k = 0; k < RUBIKS_CUBE_DIM; k++) {
                    Cube cube = rubiksCube.cubies[i][j][k];
                    for(int v = 0; v < 8; v++) 
                    {
                        cube.vertices[v] = vertex_rotate(cube.vertices[v], rotation);
                        rubiksCube.cubies[i][j][k].vertices[v] = cube.vertices[v];
                        rubiksCube.cubies[i][j][k].vertices[v] = cube.vertices[v];
                        rubiksCube.cubies[i][j][k] = updateCubeEdgesAndFaces(rubiksCube.cubies[i][j][k]);
                    }
                }
            }
        }
    }
    else if(axis == 1) {
        Matrix rotation = rotation_y(turn_spd);
        for(int i = 0; i < RUBIKS_CUBE_DIM; i++) {
            for(int j = layer; j < (layer + 1); j++) {
                for(int k = 0; k < RUBIKS_CUBE_DIM; k++) {
                    Cube cube = rubiksCube.cubies[i][j][k];
                    for(int v = 0; v < 8; v++) 
                    {
                        cube.vertices[v] = vertex_rotate(cube.vertices[v], rotation);
                        rubiksCube.cubies[i][j][k].vertices[v] = cube.vertices[v];
                        rubiksCube.cubies[i][j][k].vertices[v] = cube.vertices[v];
                        rubiksCube.cubies[i][j][k] = updateCubeEdgesAndFaces(rubiksCube.cubies[i][j][k]);
                    }
                }
            }
        }
    }
    else if(axis == 2) {
        Matrix rotation = rotation_z(turn_spd);
        for(int i = 0; i < RUBIKS_CUBE_DIM; i++) {
            for(int j = 0; j < RUBIKS_CUBE_DIM; j++) {
                for(int k = layer; k < (layer + 1); k++) {
                    Cube cube = rubiksCube.cubies[i][j][k];
                    for(int v = 0; v < 8; v++) 
                    {
                        cube.vertices[v] = vertex_rotate(cube.vertices[v], rotation);
                        rubiksCube.cubies[i][j][k].vertices[v] = cube.vertices[v];
                        rubiksCube.cubies[i][j][k].vertices[v] = cube.vertices[v];
                        rubiksCube.cubies[i][j][k] = updateCubeEdgesAndFaces(rubiksCube.cubies[i][j][k]);
                    }
                }
            }
        }
    }

    else {
        return rubiksCube;
    }

    return rubiksCube;
}








































Vertex rotateTurningVertex(Vertex v, double turn_spd, int axis) {
    if(axis == 0) {
        Matrix rotation = rotation_x(turn_spd);
        v = vertex_rotate(v, rotation);
    }
    else if(axis == 1) {
        Matrix rotation = rotation_y(turn_spd);
        v = vertex_rotate(v, rotation);
    }
    else if(axis == 2) {
        Matrix rotation = rotation_z(turn_spd);
        v = vertex_rotate(v, rotation);
    }
    return v;
}












void updateCubiesIndexing(RubiksCube *rubiksCube, int axis, int layer, int direction) 
{
    if(axis == 0) {

        Cube temp[RUBIKS_CUBE_DIM][RUBIKS_CUBE_DIM];
        Cube temp_blk[RUBIKS_CUBE_DIM][RUBIKS_CUBE_DIM]; // Temporary storage for a 2D slice

        // Copy the 2D slice into temporary storage
        for(int j = 0; j < RUBIKS_CUBE_DIM; j++) {
            for(int k = 0; k < RUBIKS_CUBE_DIM; k++) {
                temp[j][k] = rubiksCube->cubies[layer][j][k];
            }
        }


        if(direction == -1) 
        {
            for(int j = 0; j < RUBIKS_CUBE_DIM; j++) {
                for(int k = 0; k < RUBIKS_CUBE_DIM; k++) {
                    rubiksCube->cubies[layer][j][k] = temp[k][RUBIKS_CUBE_DIM-1-j];
                }
            }
        }
        else if(direction == 1) 
        {
            for(int j = 0; j < RUBIKS_CUBE_DIM; j++) {
                for(int k = 0; k < RUBIKS_CUBE_DIM; k++) {
                    rubiksCube->cubies[layer][j][k] = temp[RUBIKS_CUBE_DIM-1-k][j];
                }
            }
        }
        return;
    }

    else if(axis == 1) {

        Cube temp[RUBIKS_CUBE_DIM][RUBIKS_CUBE_DIM]; // Temporary storage for a 2D slice
        Cube temp_blk[RUBIKS_CUBE_DIM][RUBIKS_CUBE_DIM];

        // Copy the 2D slice into temporary storage
        for(int i = 0; i < RUBIKS_CUBE_DIM; i++) {
            for(int k = 0; k < RUBIKS_CUBE_DIM; k++) {
                temp[i][k] = rubiksCube->cubies[i][layer][k];
            }
        }


        if(direction == -1) 
        {
            for(int i = 0; i < RUBIKS_CUBE_DIM; i++) {
                for(int k = 0; k < RUBIKS_CUBE_DIM; k++) {
                    rubiksCube->cubies[i][layer][k] = temp[RUBIKS_CUBE_DIM-1-k][i];
                }
            }
        }
        else if(direction == 1) 
        {
            for(int i = 0; i < RUBIKS_CUBE_DIM; i++) {
                for(int k = 0; k < RUBIKS_CUBE_DIM; k++) {
                    rubiksCube->cubies[i][layer][k] = temp[k][RUBIKS_CUBE_DIM-1-i];
                }
            }
        }
        return;
    }

    else if(axis == 2) {

        Cube temp[RUBIKS_CUBE_DIM][RUBIKS_CUBE_DIM]; // Temporary storage for a 2D slice
        Cube temp_blk[RUBIKS_CUBE_DIM][RUBIKS_CUBE_DIM];

        // Copy the 2D slice into temporary storage
        for(int i = 0; i < RUBIKS_CUBE_DIM; i++) {
            for(int j = 0; j < RUBIKS_CUBE_DIM; j++) {
                temp[i][j] = rubiksCube->cubies[i][j][layer];
            }
        }


        if(direction == -1) 
        {
            for(int i = 0; i < RUBIKS_CUBE_DIM; i++) {
                for(int j = 0; j < RUBIKS_CUBE_DIM; j++) {
                    rubiksCube->cubies[i][j][layer] = temp[j][RUBIKS_CUBE_DIM-1-i];
                }
            }
        }
        else if(direction == 1) 
        {
            for(int i = 0; i < RUBIKS_CUBE_DIM; i++) {
                for(int j = 0; j < RUBIKS_CUBE_DIM; j++) {
                    rubiksCube->cubies[i][j][layer] = temp[RUBIKS_CUBE_DIM-1-j][i];
                }
            }
        }
        return;
    }


}

























































bool isOuterFace(Vertex v_rot[3], double rubiksCubeSize, bool turn_cube, bool is_turning, bool face_in_turning_layer, int axis, int layer, RubiksCube rubiksCube, int *color)
{
    Vertex v[3] = {v_rot[0], v_rot[1], v_rot[2]};
    if(is_turning == true) {
        double turning_layer_angle = rubiksCube.layer_angle[axis][layer];
        v[0] = rotateTurningVertex(v[0], -turning_layer_angle, axis);
        v[1] = rotateTurningVertex(v[1], -turning_layer_angle, axis);
        v[2] = rotateTurningVertex(v[2], -turning_layer_angle, axis);
    }

    double bias = (double)RUBIKS_CUBE_DIM / 2 - 0.01;
    bias *= rubiksCubeSize;
    int count_x = 0;
    int count_y = 0;
    int count_z = 0;


    for(int i = 0; i < 3; i++) {
            if(v[i].x < -bias || v[i].x > bias) {
                count_x++;
            }
            if(v[i].y < -bias || v[i].y > bias) {
                count_y++;
            }
            if(v[i].z < -bias || v[i].z > bias) {
                count_z++;
            }
    }
    if(count_x == 3 || count_y == 3 || count_z == 3) {
        return true;
    }
    else if(turn_cube == true && face_in_turning_layer == true) {
        if(color != NULL) {
            *color = BLACK;
        }
        return true;
    }
    return false;
}

bool isOuterFaceWithoutBlackFaces(Vertex v_rot[3], double rubiksCubeSize, bool turn_cube, bool is_turning, bool face_in_turning_layer, int axis, int layer, RubiksCube rubiksCube)
{
    Vertex v[3] = {v_rot[0], v_rot[1], v_rot[2]};
    if(is_turning == true) {
        double turning_layer_angle = rubiksCube.layer_angle[axis][layer];
        v[0] = rotateTurningVertex(v[0], -turning_layer_angle, axis);
        v[1] = rotateTurningVertex(v[1], -turning_layer_angle, axis);
        v[2] = rotateTurningVertex(v[2], -turning_layer_angle, axis);
    }

    double bias = (double)RUBIKS_CUBE_DIM / 2 - 0.01;
    bias *= rubiksCubeSize;
    int count_x = 0;
    int count_y = 0;
    int count_z = 0;


    for(int i = 0; i < 3; i++) {
            if(v[i].x < -bias || v[i].x > bias) {
                count_x++;
            }
            if(v[i].y < -bias || v[i].y > bias) {
                count_y++;
            }
            if(v[i].z < -bias || v[i].z > bias) {
                count_z++;
            }
    }
    if(count_x == 3 || count_y == 3 || count_z == 3) {
        return true;
    }
    return false;
}



bool isOuterLine(Vertex v_rot[2], double rubiksCubeSize, bool turn_cube, bool is_turning, int axis, int layer, RubiksCube rubiksCube)
{
    Vertex v[2] = {v_rot[0], v_rot[1]};
    if(is_turning == true) {
        double turning_layer_angle = rubiksCube.layer_angle[axis][layer];
        v[0] = rotateTurningVertex(v[0], -turning_layer_angle, axis);
        v[1] = rotateTurningVertex(v[1], -turning_layer_angle, axis);
    }

    double bias = (double)RUBIKS_CUBE_DIM / 2 - 0.1;
    bias *= rubiksCubeSize;
    int count_x = 0;
    int count_y = 0;
    int count_z = 0;

    for(int i = 0; i < 2; i++) {
            if(v[i].x < -bias || v[i].x > bias) {
                count_x++;
            }
            if(v[i].y < -bias || v[i].y > bias) {
                count_y++;
            }
            if(v[i].z < -bias || v[i].z > bias) {
                count_z++;
            }
    }
    if(count_x == 2 || count_y == 2 || count_z == 2) {
        return true;
    }
    return false;
}




























































Vertex crossProduct(Vertex a, Vertex b) {
    Vertex result;
    result.x = a.y * b.z - a.z * b.y;
    result.y = a.z * b.x - a.x * b.z;
    result.z = a.x * b.y - a.y * b.x;
    return result;
}

double dotProduct(Vertex a, Vertex b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vertex findNormal(Vertex v[3]) {
    Vertex AB, AC, normal;
    // Create vectors AB and AC
    AB.x = v[1].x - v[0].x;
    AB.y = v[1].y - v[0].y;
    AB.z = v[1].z - v[0].z;
    
    AC.x = v[2].x - v[0].x;
    AC.y = v[2].y - v[0].y;
    AC.z = v[2].z - v[0].z;
    
    // Compute the cross product AB × AC
    normal = crossProduct(AC, AB);
    
    return normal;
}

Vertex findNormal2(Vertex v[3]) {
    Vertex AB, AC, normal;
    // Create vectors AB and AC
    AB.x = v[1].x - v[0].x;
    AB.y = v[1].y - v[0].y;
    AB.z = v[1].z - v[0].z;
    
    AC.x = v[2].x - v[0].x;
    AC.y = v[2].y - v[0].y;
    AC.z = v[2].z - v[0].z;
    
    // Compute the cross product AC × AB
    normal = crossProduct(AB, AC);
    
    return normal;
}

Vertex normalize(Vertex v) {
    double length = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    Vertex normalized = { v.x / length, v.y / length, v.z / length };
    return normalized;
}

Vertex calculateDir(Vertex A, Vertex B) {
    Vertex dir = {
        B.x - A.x,
        B.y - A.y,
        B.z - A.z
    };
    // Normalize the direction vector
    return normalize(dir);
}

// Function to calculate the distance between two vertices
double getDistanceBetweenTwoVertices(Vertex v1, Vertex v2) {
    double x = fabs(v1.x - v2.x);
    double y = fabs(v1.y - v2.y);
    double z = fabs(v1.z - v2.z);
    return sqrt(x * x + y * y + z * z);
}

// Function to shrink a line segment by a shrink factor
void shrinkLineSegment(Vertex *P1, Vertex *P2, double shrinkFactor) {

    Vertex temp1 = *P1;
    Vertex temp2 = *P2;

    Vertex midpoint = vertex_scale(vertex_add(*P1, *P2), 0.5);

    *P1 = vertex_subtract(*P1, midpoint);
    *P2 = vertex_subtract(*P2, midpoint);

    *P1 = vertex_scale(*P1, shrinkFactor);
    *P2 = vertex_scale(*P2, shrinkFactor);

    *P1 = vertex_add(*P1, midpoint);
    *P2 = vertex_add(*P2, midpoint);
}



// various test vals for lighting
double test1 = 2.0;
double test2 = 0.05;
double test3 = 100.0;
double test4 = 0.0;
double test5 = 1.0;
double test6 = 0.5;
double test7 = 0.5;

int scaleColor(int color, int original_color, int darkest_color, double scaleFactor) {


    Uint8 a8 = color >> 24;
    Uint8 r8 = color >> 0;
    Uint8 g8 = color >> 8;
    Uint8 b8 = color >> 16;

    double a = (double)a8;
    double r = (double)r8;
    double g = (double)g8;
    double b = (double)b8;

    Uint8 oa8 = original_color >> 24;
    Uint8 or8 = original_color >> 0;
    Uint8 og8 = original_color >> 8;
    Uint8 ob8 = original_color >> 16;

    double oa = (double)oa8;
    double or = (double)or8;
    double og = (double)og8;
    double ob = (double)ob8;
    
    Uint8 da8 = darkest_color >> 24;
    Uint8 dr8 = darkest_color >> 0;
    Uint8 dg8 = darkest_color >> 8;
    Uint8 db8 = darkest_color >> 16;

    double da = (double)da8;
    double dr = (double)dr8;
    double dg = (double)dg8;
    double db = (double)db8;


    double new_r = r*scaleFactor;
    double new_g = g*scaleFactor;
    double new_b = b*scaleFactor;

    if(new_r > or) {
        new_r = or;
    }
    if(new_r < dr) {
        new_r = dr;
    }
    if(new_g > og) {
        new_g = og;
    }
    if(new_g < dg) {
        new_g = dg;
    }
    if(new_b > ob) {
        new_b = ob;
    }
    if(new_b < db) {
        new_b = db;
    }
    
    int new = (int)new_r + ((int)new_g << 8) + ((int)new_b << 16) + ((int)a << 24);

    return new;
}

int scaleColorByLightSourceDistanceAndDirection(int color, int original_color, int darkest_color, double distance_from_light_source, Vertex normal, Vertex lightDir) {

    if(distance_from_light_source == 0.0) {
        return color;
    }
    double intensity1 = 1/(distance_from_light_source * distance_from_light_source);
    intensity1 *= test3;
    intensity1 += test4;
    intensity1 = fmax(0, intensity1);

    Uint8 a8 = color >> 24;
    Uint8 r8 = color >> 0;
    Uint8 g8 = color >> 8;
    Uint8 b8 = color >> 16;

    double a = (double)a8;
    double r = (double)r8;
    double g = (double)g8;
    double b = (double)b8;

    Uint8 oa8 = original_color >> 24;
    Uint8 or8 = original_color >> 0;
    Uint8 og8 = original_color >> 8;
    Uint8 ob8 = original_color >> 16;

    double oa = (double)oa8;
    double or = (double)or8;
    double og = (double)og8;
    double ob = (double)ob8;
    
    Uint8 da8 = darkest_color >> 24;
    Uint8 dr8 = darkest_color >> 0;
    Uint8 dg8 = darkest_color >> 8;
    Uint8 db8 = darkest_color >> 16;

    double da = (double)da8;
    double dr = (double)dr8;
    double dg = (double)dg8;
    double db = (double)db8;


    Vertex normalizedNormal = normalize(normal);
    Vertex normalizedLightDir = normalize(lightDir);

    double intensity2 = dotProduct(normalizedNormal, normalizedLightDir) + test2;

    intensity2 *= test1;
    //intensity += intensity_offset;
    intensity2 = fmax(0, intensity2);  // Clamp the intensity to a minimum of 0

    double new_r = r*intensity1*intensity2;
    double new_g = g*intensity1*intensity2;
    double new_b = b*intensity1*intensity2;

    if(new_r > or) {
        new_r = or;
    }
    if(new_r < dr) {
        new_r = dr;
    }
    if(new_g > og) {
        new_g = og;
    }
    if(new_g < dg) {
        new_g = dg;
    }
    if(new_b > ob) {
        new_b = ob;
    }
    if(new_b < db) {
        new_b = db;
    }
    
    int new = (int)new_r + ((int)new_g << 8) + ((int)new_b << 16) + ((int)a << 24);

    return new;
}

double intensityOfLightSourceDistance(int color, int original_color, int darkest_color, double distance_from_light_source) {

    if(distance_from_light_source == 0.0) {
        return 1.0;
    }
    double intensity = 1/(distance_from_light_source * distance_from_light_source);
    intensity *= test3;
    intensity += test4;
    intensity = fmax(0, intensity);

    return intensity;
}

double intensityOfLightSourceDirection(int color, int original_color, int darkest_color, Vertex normal, Vertex lightDir) {


    Vertex normalizedNormal = normalize(normal);
    Vertex normalizedLightDir = normalize(lightDir);

    double intensity = dotProduct(normalizedNormal, normalizedLightDir);

    //intensity *= test1;
    //intensity += intensity_offset;
    intensity = fmax(0, intensity);  // Clamp the intensity to a minimum of 0

    return intensity;
}


















// interpolate depth in line
double interpolateDepthInLine(double x, double y, Vertex A, Vertex B) 
{

    // distances
    double dTotal = sqrt(((B.x - A.x) * (B.x - A.x)) + ((B.y - A.y) * (B.y - A.y)));
    double dPartial = sqrt(((x - A.x) * (x - A.x)) + ((y - A.y) * (y - A.y)));

    // interpolation factor
    double t = dPartial / dTotal;

    // interpolate depth value
    return A.z + (t * (B.z - A.z));
}







double signedArea(Vertex A, Vertex B, Vertex C) {
    return 0.5f * (((B.x - A.x) * (C.y - A.y)) - ((B.y - A.y) * (C.x - A.x)));
}

void ensureCCWOrder(Vertex *A, Vertex *B, Vertex *C) {
    if (signedArea(*A, *B, *C) < 0) {
        // swap A and B
        Vertex temp = *A;
        *A = *B;
        *B = temp;
    }
}

// compute barycentric coordinates of p with respect to triangle ABC
// return 1 if p is inside the triangle, 0 otherwise
int barycentric(double x, double y, Vertex A, Vertex B, Vertex C, double *u, double *v, double *w) {

    // edges
    Vertex v0 = {B.x - A.x, B.y - A.y, 0};
    Vertex v1 = {C.x - A.x, C.y - A.y, 0};
    Vertex v2 = {x - A.x, y - A.y, 0};

    // dot products
    double d00 = dotProduct(v0, v0);
    double d01 = dotProduct(v0, v1);
    double d11 = dotProduct(v1, v1);
    double d20 = dotProduct(v2, v0);
    double d21 = dotProduct(v2, v1);
    double denom = (d00 * d11) - (d01 * d01);

    if (denom == 0.0) return 0;

    *v = ((d11 * d20) - (d01 * d21)) / denom;
    *w = ((d00 * d21) - (d01 * d20)) / denom;
    *u = (1.0 - *v) - *w;

    // check if `point` is inside the triangle
    return (*u >= 0.0) && (*v >= 0.0) && (*w >= 0.0) && (*u <= 1.0) && (*v <= 1.0) && (*w <= 1.0);
}

// interpolate depth in triangle
double interpolateDepthInTri(double x, double y, Vertex A, Vertex B, Vertex C) {
    // compute barycentric coordinates
    double u, v, w;
        
    // interpolate depth value using barycentric coordinates
    double depth_in_tri;
    if(barycentric(x, y, A, B, C, &u, &v, &w)) {
        depth_in_tri =  (u * A.z) + (v * B.z) + (w * C.z);
    }
    else {
        depth_in_tri = 110.0;
    }
    return depth_in_tri;
}

double max_double(double a, double b) {
    if(a >= b) {
        return a;
    }
    else {
        return b;
    }
}

// Function to check if a line intersects a triangle
bool rayTriangleIntersectionTest(Vertex V[2], Vertex P[3]) {
    Vertex u = vertex_subtract(P[1], P[0]);
    Vertex v = vertex_subtract(P[2], P[0]);
    Vertex n = findNormal(P);
    
    if (n.x == 0 && n.y == 0 && n.z == 0) {
        return false; // Triangle is degenerate
    }

    Vertex dir = vertex_subtract(V[1], V[0]);
    Vertex w0 = vertex_subtract(V[0], P[0]);
    double a = -dotProduct(n, w0);
    double b = dotProduct(n, dir);

    if (fabs(b) == 0) {
        return false;
    }

    double r = a / b;
    if (r <= 0.0) {
        return false;
    }

    Vertex I = {V[0].x + r * dir.x, V[0].y + r * dir.y, V[0].z + r * dir.z};
    
    // Barycentric coordinates check
    double uu = dotProduct(u, u);
    double uv = dotProduct(u, v);
    double vv = dotProduct(v, v);

    Vertex w = vertex_subtract(I, P[0]);
    double wu = dotProduct(w, u);
    double wv = dotProduct(w, v);

    double D = uv * uv - uu * vv;

    double s = (uv * wv - vv * wu) / D;
    if (s < 0.0 || s > 1.0) {
        return false;
    }

    double t = (uv * wu - uu * wv) / D;
    if (t < 0.0 || (s + t) > 1.0) {
        return false;
    }

    Vertex sub_v0_I = vertex_subtract(V[0], I);
    double intersectionDist = sqrt(dotProduct(sub_v0_I, sub_v0_I));
    double pixelAtSurfaceDist = sqrt(dotProduct(dir, dir));
    
    if(intersectionDist < pixelAtSurfaceDist) {
        return true;
    }
    
    return false;
}

bool shadowPixel(Vertex pixelPositionAtSurface, Vertex lightSourcePosition, uint32_t *objectAddress, double *intensity) {

    Vertex ray[2] = {pixelPositionAtSurface, lightSourcePosition};
    for(int i = 0; i < collision_n; i++) {

        if(collisionTris.object_address[i] == objectAddress) {
            continue;
        }

        Vertex tri[3] = {collisionTris.v[i][0], collisionTris.v[i][1], collisionTris.v[i][2]};

        if(rayTriangleIntersectionTest(ray, tri)) {
            *intensity = 0.3;
            return true;
        }
    }
    return false;
}





Vertex reflect(Vertex incident, Vertex normal) {
    normal = normalize(normal);
    // Calculate the dot product of the incident vector and the normal
    double dotProd = dotProduct(incident, normal);

    // Reflect the incident vector off the surface
    Vertex reflection;
    reflection.x = incident.x - 2 * dotProd * normal.x;
    reflection.y = incident.y - 2 * dotProd * normal.y;
    reflection.z = incident.z - 2 * dotProd * normal.z;

    return reflection;
}
















typedef struct {
    Vertex start;
    Vertex end;
} Edge;

void swapEdges(Edge* v1, Edge* v2) {
    Edge temp = *v1;
    *v1 = *v2;
    *v2 = temp;
}

double computeXIntersection(Edge edge, double y) {
    if ((int)edge.start.y == (int)edge.end.y) {
        return edge.end.x;
    } 
    else {
        double t = (double)(y - (int)edge.start.y) / (double)((int)edge.end.y - (int)edge.start.y);
        return (edge.start.x + t * (edge.end.x - edge.start.x));
    }
}


typedef struct {
    Edge edges[3];
    int count;
} EdgeList;

void updateAEL(int y, EdgeList* AEL, Edge edges[3]) {
    int i = 0;
    AEL->count = 0;
    for (int j = 0; j < 3; j++) {
        if (((int)edges[j].start.y <= y && (int)edges[j].end.y > y) || ((int)edges[j].start.y > y && (int)edges[j].end.y <= y)) {
            AEL->edges[AEL->count] = edges[j];
            AEL->count++;
        }
        if(AEL->count >= 2) {
            break;
        }

    }
    return;
}


#define PIXELS ((WINDOW_HEIGHT) * (WINDOW_WIDTH))

typedef struct {
    int color;
    double depth;
} Pixel;




//const double denom = 20.0;
double far_distance = 50;
double coords_scalar = 0.8;
uint32_t *floorGridObjPtr;


void getTriPixels(Pixel **pixels, Vertex p[3], Vertex v[3], Vertex surfaceNormal, Vertex surfacePosition, int color, Vertex lightSourcePosition, Vertex camera, Matrix invRotX, Matrix invRotY, uint32_t *objectAddress) 
{
    Vertex surfaceNormalNormalized = normalize(findNormal(v));
    
    
    Vertex tempv = center_of_tri(v);

    v[0] = vertex_subtract(v[0], tempv);
    v[1] = vertex_subtract(v[1], tempv);
    v[2] = vertex_subtract(v[2], tempv);

    v[0] = vertex_scale(v[0], 2.0);
    v[1] = vertex_scale(v[1], 2.0);
    v[2] = vertex_scale(v[2], 2.0);

    v[0] = vertex_add(v[0], tempv);
    v[1] = vertex_add(v[1], tempv);
    v[2] = vertex_add(v[2], tempv);

    
    // y min
    double minY = p[0].y;
    for(int i = 1; i < 3; i++) {
        if(p[i].y < minY) {
            minY = p[i].y;
        }
    }
    // y max
    double maxY = p[0].y;
    for(int i = 1; i < 3; i++) {
        if(p[i].y > maxY) {
            maxY = p[i].y;
        }
    }
    
    Edge edges[3] = {
        {.start = p[0], .end = p[1]},
        {.start = p[1], .end = p[2]},
        {.start = p[2], .end = p[0]},
    };

    EdgeList AEL = {.count = 0};

    for (int y = (int)minY; y <= (int)maxY; y++) {
        if ((int)y < 1 || (int)y >= WINDOW_HEIGHT-2) {
            continue;
        }
        updateAEL(y, &AEL, edges);

        double xValues[3];
        int xCount = 0;

        for (int i = 0; i < AEL.count; i++) {
            xValues[xCount] = computeXIntersection(AEL.edges[i], y);
            xCount++;
        }

        if (xValues[0] > xValues[1]) {
            double temp = xValues[0];
            xValues[0] = xValues[1];
            xValues[1] = temp;
        }

        double startX = xValues[0];
        double endX = xValues[1];

        for (double x = startX; x <= endX; x++) {

            if ((int)x < 1 || (int)x >= WINDOW_WIDTH-2) {
                continue;
            }
            
            double depth;
            double px;
            double py;
            inv_view_transform_vertex_f(x, y, &px, &py);
            depth = interpolateDepthInTri(px, py, v[0], v[1], v[2]);
            if (depth < pixels[(int)x][(int)y].depth) {
                
                Vertex inv_view = {px, py, depth};
                Vertex a = vertex_project(inv_view, deprojection);
                Vertex deprojected = vertex_rotate(a, invRotX);
                deprojected = vertex_rotate(deprojected, invRotY);
                deprojected.x += camera.x;
                deprojected.y += camera.y;
                deprojected.z += camera.z;

                double distance_from_player = getDistanceBetweenTwoVertices(deprojected, camera);
                if(distance_from_player > far_distance) {
                    continue;
                }

                pixels[(int)x][(int)y].depth = depth;
                
                double lightDist = getDistanceBetweenTwoVertices(deprojected, lightSourcePosition);
                double intensityOfLightDist = 2*test3 / (lightDist * lightDist);

                double intensityOfCamDist = 10.0 / (distance_from_player);

                Vertex lightDir = normalize(calculateDir(deprojected, lightSourcePosition)); // from `a` to light source
                double cosTheta = dotProduct(normalize(surfaceNormal), lightDir);
                double intensityOfLightDir = fmax(0.0, cosTheta) + test2;
                

                // stuff for reflective surfaces
                /*
                Vertex viewDir = calculateDir(deprojected, camera); // from `a` to pixel 
                double shininess = test5;
                double diffuse = test6;
                double specular = test7;
                Vertex revLightDir = vertex_scale(lightDir, -1.0);
                Vertex reflectDir = reflect(revLightDir, surfaceNormalNormalized);
                double spec = pow(fmax(dotProduct(viewDir, reflectDir), 0.0), shininess);
                double intensityOfReflection = (diffuse * cosTheta) + (specular * spec);
                */
                
                
                int darkest_color = 0xFF000000;
                int colorA = scaleColor(color, color, darkest_color, (intensityOfLightDist*intensityOfCamDist*intensityOfLightDir));
                if(objectAddress != NULL) {
                    
                    double intensityOfShadow = 1.0;
                    bool isShadowPixel;
                    if(drawShadows == false) {
                        isShadowPixel = false;
                    }
                    else if(drawShadows = true) {
                        isShadowPixel = shadowPixel(deprojected, lightSourcePosition, objectAddress, &intensityOfShadow);
                    }
                    if(isShadowPixel == true) {
                        int colorB = scaleColor(colorA, color, darkest_color, intensityOfShadow);
                        pixels[(int)x][(int)y].color = colorB;
                    }

                    else {
                        pixels[(int)x][(int)y].color = colorA;
                    }
                }

                else if(objectAddress == NULL) {
                    pixels[(int)x][(int)y].color = colorA;
                }
                
            }        
        }
    }
}

















void thickenPixel(Pixel **pixels, int x, int y, double depth, int color) {
    double depths[10] = {0.83, 0.81, 0.75, 0.7, 0.65, 0.6, 0.5, 0.4, 0.3, 0.2};
    //for(int i = 0; i < 10; i++) {
        //depths[i] = 1.0;
    //}
    if (depth < pixels[x][y].depth) 
    {
        // Center
        pixels[x][y].color = color;
        //pixels[x][y].depth = depth;

        if(depth < depths[0]) {
            if(depth < pixels[x][y-1].depth) {
                pixels[x][y-1].color = color;
                //pixels[x][y].depth = depth;
            }
        }

        if(depth < depths[1]) {
            if(depth < pixels[x-1][y].depth) {
                pixels[x-1][y].color = color;
                //pixels[x-1][y].depth = depth;
            }
        }

        if(depth < depths[2]) {
            if(depth < pixels[x+1][y].depth) {
                pixels[x+1][y].color = color;
                //pixels[x+1][y].depth = depth;
            }
            if(depth < pixels[x][y+1].depth) {
                pixels[x][y+1].color = color;
                //pixels[x][y+1].depth = depth;
            }
        }

        if(depth < depths[3]) {
            if(depth < pixels[x-1][y-1].depth) {
                pixels[x-1][y-1].color = color;
                //pixels[x-1][y-1].depth = depth;
            }
            if(depth < pixels[x+1][y+1].depth) {
                pixels[x+1][y+1].color = color;
                //pixels[x+1][y+1].depth = depth;
            }
        }

        if(depth < depths[4]) {
            if(depth < pixels[x+1][y-1].depth) {
                pixels[x+1][y-1].color = color;
                //pixels[x+1][y-1].depth = depth;
            }
            if(depth < pixels[x-1][y+1].depth) {
                pixels[x-1][y+1].color = color;
                //pixels[x-1][y+1].depth = depth;
            }
        }

        if(depth < depths[5]) {
            if(depth < pixels[x-2][y].depth) {
                pixels[x-2][y].color = color;
                //pixels[x-2][y].depth = depth;
            }
            if(depth < pixels[x][y-2].depth) {
                pixels[x][y-2].color = color;
                //pixels[x][y-2].depth = depth;
            }
        }
        if(depth < depths[6]) {
            if(depth < pixels[x+2][y].depth) {
                pixels[x+2][y].color = color;
                //pixels[x+2][y].depth = depth;
            }
            if(depth < pixels[x][y+2].depth) {
                pixels[x][y+2].color = color;
                //pixels[x][y+2].depth = depth;
            }
        }
        if(depth < depths[7]) {
            if(depth < pixels[x-2][y-1].depth) { 
                pixels[x-2][y-1].color = color;
                //pixels[x-2][y-1].depth = depth;
            }
            if(depth < pixels[x+2][y-1].depth) {
                pixels[x+2][y-1].color = color;
                //pixels[x+2][y-1].depth = depth;
            }
            if(depth < pixels[x-2][y+1].depth) {
                pixels[x-2][y+1].color = color;
                //pixels[x-2][y+1].depth = depth;
            }
            if(depth < pixels[x+2][y+1].depth) {
                pixels[x+2][y+1].color = color;
                //pixels[x+2][y+1].depth = depth;
            }
        }
        if(depth < depths[8]) {
            if(depth < pixels[x-1][y-2].depth) {
                pixels[x-1][y-2].color = color;
                //pixels[x-1][y-2].depth = depth;
            }
            if(depth < pixels[x+1][y-2].depth) {
                pixels[x+1][y-2].color = color;
                //pixels[x+1][y-2].depth = depth;
            }
            if(depth < pixels[x-1][y+2].depth) {
                pixels[x-1][y+2].color = color;
                //pixels[x-1][y+2].depth = depth;
            }
            if(depth < pixels[x+1][y+2].depth) {
                pixels[x+1][y+2].color = color;
                //pixels[x+1][y+2].depth = depth;
            }
        }
        if(depth < depths[9]) {
            if(depth < pixels[x-2][y-2].depth) {
                pixels[x-2][y-2].color = color;
                //pixels[x-2][y-2].depth = depth;
            }
            if(depth < pixels[x+2][y-2].depth) {
                pixels[x+2][y-2].color = color;
                //pixels[x+2][y-2].depth = depth;
            }
            if(depth < pixels[x-2][y+2].depth) {
                pixels[x-2][y+2].color = color;
                //pixels[x-2][y+2].depth = depth;
            }
            if(depth < pixels[x+2][y+2].depth) {
                pixels[x+2][y+2].color = color;
                //pixels[x+2][y+2].depth = depth;
            }
        }
    }
}

void getLinePixels(Pixel **pixels, Vertex p[2], Vertex v[2], int color, Vertex camera, Matrix invRotX, Matrix invRotY, uint32_t *objectAddress)
{
    double x, y, dx, dy, dx1, dy1, px, py, xe, ye, i;
    dx = p[1].x - p[0].x; 
    dy = p[1].y - p[0].y;

    dx1 = fabs(dx);
    dy1 = fabs(dy);

    px = 2 * dy1 - dx1; 
    py = 2 * dx1 - dy1;

    if(dy1 <= dx1) 
    {
        if(dx >= 0) 
        {
            x = p[0].x; 
            y = p[0].y; 
            xe = p[1].x;
        }
        else 
        {
            x = p[1].x; 
            y = p[1].y; 
            xe = p[0].x;
        }

        double fx;
        double fy;
        inv_view_transform_vertex_f(x, y, &fx, &fy);
        double depth = interpolateDepthInLine(fx, fy, v[0], v[1]);

        Vertex inv_view = {fx, fy, depth};
        Vertex a = vertex_project(inv_view, deprojection);
        Vertex deprojected = vertex_rotate(a, invRotX);
        deprojected = vertex_rotate(deprojected, invRotY);
        deprojected.x += camera.x;
        deprojected.y += camera.y;
        deprojected.z += camera.z;

        double distance_from_player = getDistanceBetweenTwoVertices(deprojected, camera);

        if(!((int)x < 2 || (int)x >= WINDOW_WIDTH-2 || (int)y < 2 || (int)y >= WINDOW_HEIGHT-2)) {
            if(distance_from_player < far_distance) {
                thickenPixel(pixels, (int)x, (int)y, depth, color);
            }
        }

        for(i = 0; x < xe; i++) 
        {
            x = x + 1;
            if(px < 0) 
            {
                px = px + 2 * dy1;
            }
            else 
            {
                if((dx < 0 && dy < 0) || (dx > 0 && dy > 0)) 
                {
                    y = y + 1;
                }
                else 
                {
                    y = y - 1;
                }
                px = px + 2 * (dy1 - dx1);
            }

        double fx;
        double fy;
        inv_view_transform_vertex_f(x, y, &fx, &fy);
        double depth = interpolateDepthInLine(fx, fy, v[0], v[1]);

        Vertex inv_view = {fx, fy, depth};
        Vertex a = vertex_project(inv_view, deprojection);
        Vertex deprojected = vertex_rotate(a, invRotX);
        deprojected = vertex_rotate(deprojected, invRotY);
        deprojected.x += camera.x;
        deprojected.y += camera.y;
        deprojected.z += camera.z;

        double distance_from_player = getDistanceBetweenTwoVertices(deprojected, camera);

        if(!((int)x < 2 || (int)x >= WINDOW_WIDTH-2 || (int)y < 2 || (int)y >= WINDOW_HEIGHT-2)) {
            if(distance_from_player < far_distance) {
                thickenPixel(pixels, (int)x, (int)y, depth, color);
            }
        }
        }
    }
    else 
    {
        if(dy >= 0) 
        {
            x = p[0].x; 
            y = p[0].y; 
            ye = p[1].y;
        }
        else 
        {
            x = p[1].x; 
            y = p[1].y; 
            ye = p[0].y;
        }

        double fx;
        double fy;
        inv_view_transform_vertex_f(x, y, &fx, &fy);
        double depth = interpolateDepthInLine(fx, fy, v[0], v[1]);

        Vertex inv_view = {fx, fy, depth};
        Vertex a = vertex_project(inv_view, deprojection);
        Vertex deprojected = vertex_rotate(a, invRotX);
        deprojected = vertex_rotate(deprojected, invRotY);
        deprojected.x += camera.x;
        deprojected.y += camera.y;
        deprojected.z += camera.z;

        double distance_from_player = getDistanceBetweenTwoVertices(deprojected, camera);

        if(!((int)x < 2 || (int)x >= WINDOW_WIDTH-2 || (int)y < 2 || (int)y >= WINDOW_HEIGHT-2)) {
            if(distance_from_player < far_distance) {
                thickenPixel(pixels, (int)x, (int)y, depth, color);
            }
        }

        for(i = 0.0; y < ye; i++) 
        {
            y = y + 1;
            if(py <= 0) 
            {
                py = py + 2 * dx1;
            }
            else 
            {
                if((dx < 0 && dy < 0) || (dx > 0 && dy > 0)) 
                {
                    x = x + 1;
                }
                else 
                {
                    x = x - 1;
                }
                py = py + 2 * (dx1 - dy1);
            }

                double fx;
        double fy;
        inv_view_transform_vertex_f(x, y, &fx, &fy);
        double depth = interpolateDepthInLine(fx, fy, v[0], v[1]);

        Vertex inv_view = {fx, fy, depth};
        Vertex a = vertex_project(inv_view, deprojection);
        Vertex deprojected = vertex_rotate(a, invRotX);
        deprojected = vertex_rotate(deprojected, invRotY);
        deprojected.x += camera.x;
        deprojected.y += camera.y;
        deprojected.z += camera.z;

        double distance_from_player = getDistanceBetweenTwoVertices(deprojected, camera);

        if(!((int)x < 2 || (int)x >= WINDOW_WIDTH-2 || (int)y < 2 || (int)y >= WINDOW_HEIGHT-2)) {
            if(distance_from_player < far_distance) {
                thickenPixel(pixels, (int)x, (int)y, depth, color);
            }
        }
        }
    }
}















































































int isCounterClockwiseTri(Vertex v[3]) {
    int crossProduct = (v[1].x - v[0].x) * (v[2].y - v[0].y) - (v[1].y - v[0].y) * (v[2].x - v[0].x);

    if (crossProduct > 0) {
        return 1;  // Points are in counterclockwise order
    } else if (crossProduct < 0) {
        return -1; // Points are in clockwise order
    } else {
        return 0;  // Points are collinear
    }
}














typedef struct {
    Vertex vertices[4];
    Vertex edges[4][2];
    Face faces[2];
} Square;

Square createSquare(Vertex v, double s, int color) {
    
    Square square = {
        .vertices = {
            {-0.5 * s + v.x, 0.0 + v.y, -0.5 * s + v.z},
            {0.5 * s + v.x, 0.0 + v.y, -0.5 * s + v.z},
            {0.5 * s + v.x, 0.0 + v.y, 0.5 * s + v.z},
            {-0.5 * s + v.x, 0.0 + v.y, 0.5 * s + v.z},
        },
        .edges = {
            {square.vertices[0], square.vertices[1]}, // 0, 1
            {square.vertices[1], square.vertices[2]}, // 1, 2
            {square.vertices[2], square.vertices[3]}, // 2, 3
            {square.vertices[3], square.vertices[0]}, // 3, 0
        },
        .faces[0] = {   
            .tri = {square.vertices[0], square.vertices[1], square.vertices[2]},
            .color = color
        },
        .faces[1] = {   
            .tri = {square.vertices[2], square.vertices[3], square.vertices[0]},
            .color = color
        },
    };
    return square;
}

int FLOOR_GRID_DIM = 11;

typedef struct {
    Square **squares;
} FloorGrid;

FloorGrid createFloorGrid(FloorGrid floorGrid, double square_scale, int color) {

    floorGrid.squares = malloc((FLOOR_GRID_DIM) * sizeof(Square*));
    for(int i = 0; i < FLOOR_GRID_DIM; i++) {
        floorGrid.squares[i] = malloc((FLOOR_GRID_DIM) * sizeof(Square));
    }

    double offset = ((double)FLOOR_GRID_DIM / 2) - 0.5;
    Vertex v;
    for(int i = 0; i < FLOOR_GRID_DIM; i++) {
        for(int j = 0; j < FLOOR_GRID_DIM; j++) {
            v.x = ((double)i - offset) * square_scale;
            v.z = ((double)j - offset) * square_scale;
            v.y = -.05;//-2.77;
            floorGrid.squares[i][j] = createSquare(v, square_scale, color);
        }
    }
    return floorGrid;
}






bool isFaceInTurningLayer(Vertex v_scaled[3], int axis, int layer, Vertex rubiksCubePosition, double rubiksCubeSize) {
    Vertex v[3] = {v_scaled[0], v_scaled[1], v_scaled[2]};
    
    for(int i = 0; i < 3; i++) {
        v[i].x /= rubiksCubeSize;
        v[i].y /= rubiksCubeSize;
        v[i].z /= rubiksCubeSize;
    }
    double slices[3];
    slices[0] = ((layer + 0.5) - ((double)RUBIKS_CUBE_DIM/2)) - 1;
    slices[1] = ((layer + 0.5) - ((double)RUBIKS_CUBE_DIM/2));
    slices[2] = ((layer + 0.5) - ((double)RUBIKS_CUBE_DIM/2)) + 1;

    switch(axis) {
        case 0:
            if(v[0].x > slices[1] && v[1].x > slices[1] && v[2].x > slices[1] && v[0].x < slices[2] && v[1].x < slices[2] && v[2].x < slices[2] || v[0].x < slices[1] && v[1].x < slices[1] && v[2].x < slices[1] && v[0].x > slices[0] && v[1].x > slices[0] && v[2].x > slices[0]) {
                return true;
            }
            break;
        case 1:
            if(v[0].y > slices[1] && v[1].y > slices[1] && v[2].y > slices[1] && v[0].y < slices[2] && v[1].y < slices[2] && v[2].y < slices[2] || v[0].y < slices[1] && v[1].y < slices[1] && v[2].y < slices[1] && v[0].y > slices[0] && v[1].y > slices[0] && v[2].y > slices[0]) {
                return true;
            }
            break;
        case 2:
            if(v[0].z > slices[1] && v[1].z > slices[1] && v[2].z > slices[1] && v[0].z < slices[2] && v[1].z < slices[2] && v[2].z < slices[2] || v[0].z < slices[1] && v[1].z < slices[1] && v[2].z < slices[1] && v[0].z > slices[0] && v[1].z > slices[0] && v[2].z > slices[0]) {
                return true;
            }
            break;
    }
    return false;
}



















/*int clip_tri(Vertex v[3], Vertex tri[2][3]) {

    Vertex center = center_of_polygon(v, 3);
    double x = center.x / center.z;
    double y = center.y / center.z;

    // start of near clipping
    Vertex temp0[3] = {v[0], v[1], v[2]};
    Vertex temp1[3] = {v[0], v[1], v[2]};
    Vertex temp2[3] = {v[0], v[1], v[2]};

    int near_clip[3];
    near_clip[0] = clip_vertex_near(&temp0[0], &temp0[1]);
    near_clip[1] = clip_vertex_near(&temp1[1], &temp1[2]);
    near_clip[2] = clip_vertex_near(&temp2[2], &temp2[0]);

    // if all 3 vertices are in front of the near plane, then skip tri
    if(near_clip[0] == 0 && near_clip[1] == 0 && near_clip[2] == 0) {
        return 0;
    }

    int tri_count = 0;
    bool needs_near_clipping = true;
    bool check_far = false;

    // if all 3 vertices are behind the near plane, no clipping needed
    if(near_clip[0] == 3 && near_clip[1] == 3 && near_clip[2] == 3) {
        check_far = true;
        tri_count = 1;
        needs_near_clipping = false;
    }

    // otherwise, create new triangle(s) if necessary
    if(needs_near_clipping = true) {
        if(near_clip[0] == 1 && near_clip[1] == 3 && near_clip[2] == 2) { // 2
            tri[tri_count][0] = temp0[0];
            tri[tri_count][1] = temp0[1];
            tri[tri_count][2] = temp0[2];
            tri_count++;
            tri[tri_count][0] = temp0[2];
            tri[tri_count][1] = temp2[0];
            tri[tri_count][2] = temp0[0];
            tri_count++;                                
        }
        else if(near_clip[0] == 0 && near_clip[1] == 1 && near_clip[2] == 2) { // 1
            tri[tri_count][0] = temp2[0];
            tri[tri_count][1] = temp1[1];
            tri[tri_count][2] = temp1[2];
            tri_count++;
        }
        else if(near_clip[0] == 1 && near_clip[1] == 2 && near_clip[2] == 0) { // 1
            tri[tri_count][0] = temp0[0];
            tri[tri_count][1] = temp0[1];
            tri[tri_count][2] = temp1[2];
            tri_count++;
        }
        else if(near_clip[0] == 2 && near_clip[1] == 1 && near_clip[2] == 3) { // 2 
            tri[tri_count][0] = temp0[0];
            tri[tri_count][1] = temp0[1];
            tri[tri_count][2] = temp2[2];
            tri_count++;
            tri[tri_count][0] = temp1[1];
            tri[tri_count][1] = temp1[2];
            tri[tri_count][2] = temp0[1];
            tri_count++;
        }
        else if(near_clip[0] == 2 && near_clip[1] == 0 && near_clip[2] == 1) { // 1
            tri[tri_count][0] = temp0[0];
            tri[tri_count][1] = temp0[1];
            tri[tri_count][2] = temp2[2];
            tri_count++;
        }
        else if(near_clip[0] == 3 && near_clip[1] == 2 && near_clip[2] == 1) { // 2
            tri[tri_count][0] = temp0[0];
            tri[tri_count][1] = temp0[1];
            tri[tri_count][2] = temp1[2];
            tri_count++;
            tri[tri_count][0] = temp1[2];
            tri[tri_count][1] = temp2[2];
            tri[tri_count][2] = temp0[0];
            tri_count++;
        }
        else {
            tri_count++;
        }
    }

    // start of far clipping
    if(check_far == true) {

        Vertex temp3[3] = {v[0], v[1], v[2]};
        Vertex temp4[3] = {v[0], v[1], v[2]};
        Vertex temp5[3] = {v[0], v[1], v[2]};

        int far_clip[3];
        far_clip[0] = clip_vertex_far(&temp3[0], &temp3[1]);
        far_clip[1] = clip_vertex_far(&temp4[1], &temp4[2]);
        far_clip[2] = clip_vertex_far(&temp5[2], &temp5[0]);

        bool needs_far_clipping = true;
        // if all 3 vertices are behind of the far plane, then skip tri
        if(far_clip[0] == 0 && far_clip[1] == 0 && far_clip[2] == 0) {
            return 0;
        }
        // if all 3 vertices are in front of the far plane, then no far clipping is needed
        else if(far_clip[0] == 3 && far_clip[1] == 3 && far_clip[2] == 3) {
            tri_count = 1;
            needs_far_clipping = false;
        }


        if(needs_far_clipping == true) {
            tri_count = 0;
            if(far_clip[0] == 1 && far_clip[1] == 3 && far_clip[2] == 2) { // 2
                tri[tri_count][0] = temp3[0];
                tri[tri_count][1] = temp3[1];
                tri[tri_count][2] = temp3[2];
                tri_count++;
                tri[tri_count][0] = temp3[2];
                tri[tri_count][1] = temp5[0];
                tri[tri_count][2] = temp3[0];
                tri_count++;                                
            }
            else if(far_clip[0] == 0 && far_clip[1] == 1 && far_clip[2] == 2) { // 1
                tri[tri_count][0] = temp5[0];
                tri[tri_count][1] = temp4[1];
                tri[tri_count][2] = temp4[2];
                tri_count++;
            }
            else if(far_clip[0] == 1 && far_clip[1] == 2 && far_clip[2] == 0) { // 1
                tri[tri_count][0] = temp3[0];
                tri[tri_count][1] = temp3[1];
                tri[tri_count][2] = temp4[2];
                tri_count++;
            }
            else if(far_clip[0] == 2 && far_clip[1] == 1 && far_clip[2] == 3) { // 2 
                tri[tri_count][0] = temp3[0];
                tri[tri_count][1] = temp3[1];
                tri[tri_count][2] = temp5[2];
                tri_count++;
                tri[tri_count][0] = temp4[1];
                tri[tri_count][1] = temp4[2];
                tri[tri_count][2] = temp3[1];
                tri_count++;
            }
            else if(far_clip[0] == 2 && far_clip[1] == 0 && far_clip[2] == 1) { // 1
                tri[tri_count][0] = temp3[0];
                tri[tri_count][1] = temp3[1];
                tri[tri_count][2] = temp5[2];
                tri_count++;
            }
            else if(far_clip[0] == 3 && far_clip[1] == 2 && far_clip[2] == 1) { // 2
                tri[tri_count][0] = temp3[0];
                tri[tri_count][1] = temp3[1];
                tri[tri_count][2] = temp4[2];
                tri_count++;
                tri[tri_count][0] = temp4[2];
                tri[tri_count][1] = temp5[2];
                tri[tri_count][2] = temp3[0];
                tri_count++;
            }
            else {
                tri_count++;
            }
        }
    }
    return tri_count;
}*/

int clip_tri(Vertex v[3], Vertex tri[2][3]) {

// return value
    int tri_count = 0;

// get clipped edges
    Vertex temp0[3] = {v[0], v[1], v[2]};
    Vertex temp1[3] = {v[0], v[1], v[2]};
    Vertex temp2[3] = {v[0], v[1], v[2]};

    int near_clip[3];
    near_clip[0] = clip_vertex_near(&temp0[0], &temp0[1]);
    near_clip[1] = clip_vertex_near(&temp1[1], &temp1[2]);
    near_clip[2] = clip_vertex_near(&temp2[2], &temp2[0]);

    Vertex temp3[3] = {v[0], v[1], v[2]};
    Vertex temp4[3] = {v[0], v[1], v[2]};
    Vertex temp5[3] = {v[0], v[1], v[2]};

    int far_clip[3];
    far_clip[0] = clip_vertex_far(&temp3[0], &temp3[1]);
    far_clip[1] = clip_vertex_far(&temp4[1], &temp4[2]);
    far_clip[2] = clip_vertex_far(&temp5[2], &temp5[0]);

// check if all 3 vertices are outside or inside of near or far plane
    bool outside_of_near = near_clip[0] == 0 && near_clip[1] == 0 && near_clip[2] == 0;
    bool outside_of_far = far_clip[0] == 0 && far_clip[1] == 0 && far_clip[2] == 0;
    // use xor (if both are true then the triangle is probably visible)
    if(outside_of_near ^ outside_of_far) {
        return 0; // draw 0 tris
    }

    bool inside_of_near = near_clip[0] == 3 && near_clip[1] == 3 && near_clip[2] == 3;
    bool inside_of_far = far_clip[0] == 3 && far_clip[1] == 3 && far_clip[2] == 3;
    if(inside_of_near && inside_of_far) {
        return 1; // draw 1 tri
    }

// otherwise, create new triangle(s) if necessary
    if(inside_of_near == false) {
        if(near_clip[0] == 1 && near_clip[1] == 3 && near_clip[2] == 2) { // 2
            tri[tri_count][0] = temp0[0];
            tri[tri_count][1] = temp0[1];
            tri[tri_count][2] = temp0[2];
            tri_count++;
            tri[tri_count][0] = temp0[2];
            tri[tri_count][1] = temp2[0];
            tri[tri_count][2] = temp0[0];
            tri_count++;                                
        }
        else if(near_clip[0] == 0 && near_clip[1] == 1 && near_clip[2] == 2) { // 1
            tri[tri_count][0] = temp2[0];
            tri[tri_count][1] = temp1[1];
            tri[tri_count][2] = temp1[2];
            tri_count++;
        }
        else if(near_clip[0] == 1 && near_clip[1] == 2 && near_clip[2] == 0) { // 1
            tri[tri_count][0] = temp0[0];
            tri[tri_count][1] = temp0[1];
            tri[tri_count][2] = temp1[2];
            tri_count++;
        }
        else if(near_clip[0] == 2 && near_clip[1] == 1 && near_clip[2] == 3) { // 2 
            tri[tri_count][0] = temp0[0];
            tri[tri_count][1] = temp0[1];
            tri[tri_count][2] = temp2[2];
            tri_count++;
            tri[tri_count][0] = temp1[1];
            tri[tri_count][1] = temp1[2];
            tri[tri_count][2] = temp0[1];
            tri_count++;
        }
        else if(near_clip[0] == 2 && near_clip[1] == 0 && near_clip[2] == 1) { // 1
            tri[tri_count][0] = temp0[0];
            tri[tri_count][1] = temp0[1];
            tri[tri_count][2] = temp2[2];
            tri_count++;
        }
        else if(near_clip[0] == 3 && near_clip[1] == 2 && near_clip[2] == 1) { // 2
            tri[tri_count][0] = temp0[0];
            tri[tri_count][1] = temp0[1];
            tri[tri_count][2] = temp1[2];
            tri_count++;
            tri[tri_count][0] = temp1[2];
            tri[tri_count][1] = temp2[2];
            tri[tri_count][2] = temp0[0];
            tri_count++;
        }
        else {
            tri_count++;
        }
    }

    

    // far clip only
    else if(inside_of_far == false && inside_of_near == true) {
        if(far_clip[0] == 1 && far_clip[1] == 3 && far_clip[2] == 2) { // 2
            tri[tri_count][0] = temp3[0];
            tri[tri_count][1] = temp3[1];
            tri[tri_count][2] = temp3[2];
            tri_count++;
            tri[tri_count][0] = temp3[2];
            tri[tri_count][1] = temp5[0];
            tri[tri_count][2] = temp3[0];
            tri_count++;                                
        }
        else if(far_clip[0] == 0 && far_clip[1] == 1 && far_clip[2] == 2) { // 1
            tri[tri_count][0] = temp5[0];
            tri[tri_count][1] = temp4[1];
            tri[tri_count][2] = temp4[2];
            tri_count++;
        }
        else if(far_clip[0] == 1 && far_clip[1] == 2 && far_clip[2] == 0) { // 1
            tri[tri_count][0] = temp3[0];
            tri[tri_count][1] = temp3[1];
            tri[tri_count][2] = temp4[2];
            tri_count++;
        }
        else if(far_clip[0] == 2 && far_clip[1] == 1 && far_clip[2] == 3) { // 2 
            tri[tri_count][0] = temp3[0];
            tri[tri_count][1] = temp3[1];
            tri[tri_count][2] = temp5[2];
            tri_count++;
            tri[tri_count][0] = temp4[1];
            tri[tri_count][1] = temp4[2];
            tri[tri_count][2] = temp3[1];
            tri_count++;
        }
        else if(far_clip[0] == 2 && far_clip[1] == 0 && far_clip[2] == 1) { // 1
            tri[tri_count][0] = temp3[0];
            tri[tri_count][1] = temp3[1];
            tri[tri_count][2] = temp5[2];
            tri_count++;
        }
        else if(far_clip[0] == 3 && far_clip[1] == 2 && far_clip[2] == 1) { // 2
            tri[tri_count][0] = temp3[0];
            tri[tri_count][1] = temp3[1];
            tri[tri_count][2] = temp4[2];
            tri_count++;
            tri[tri_count][0] = temp4[2];
            tri[tri_count][1] = temp5[2];
            tri[tri_count][2] = temp3[0];
            tri_count++;
        }
        else {
            tri_count++;
        }
    }

    // both planes
    else if(inside_of_near == false && inside_of_far == false) {
        int n = tri_count;
        for(int i = 0; i < n; i++) {

            Vertex temp3b[3] = {tri[i][0], tri[i][1], tri[i][2]};
            Vertex temp4b[3] = {tri[i][0], tri[i][1], tri[i][2]};
            Vertex temp5b[3] = {tri[i][0], tri[i][1], tri[i][2]};

            far_clip[0] = clip_vertex_far(&temp3b[0], &temp3b[1]);
            far_clip[1] = clip_vertex_far(&temp4b[1], &temp4b[2]);
            far_clip[2] = clip_vertex_far(&temp5b[2], &temp5b[0]);

            if(far_clip[0] == 3 && far_clip[1] == 3 && far_clip[2] == 3) { // skip
                continue;
            }

            if(far_clip[0] == 1 && far_clip[1] == 3 && far_clip[2] == 2) { // split
                tri[i][0] = temp3b[0];
                tri[i][1] = temp3b[1];
                tri[i][2] = temp3b[2];

                tri[tri_count][0] = temp3b[2];
                tri[tri_count][1] = temp5b[0];
                tri[tri_count][2] = temp3b[0];
                tri_count++;                                
            }
            else if(far_clip[0] == 0 && far_clip[1] == 1 && far_clip[2] == 2) { // dont split
                tri[i][0] = temp5b[0];
                tri[i][1] = temp4b[1];
                tri[i][2] = temp4b[2];

            }
            else if(far_clip[0] == 1 && far_clip[1] == 2 && far_clip[2] == 0) { // dont split
                tri[i][0] = temp3b[0];
                tri[i][1] = temp3b[1];
                tri[i][2] = temp4b[2];

            }
            else if(far_clip[0] == 2 && far_clip[1] == 1 && far_clip[2] == 3) { // split
                tri[i][0] = temp3b[0];
                tri[i][1] = temp3b[1];
                tri[i][2] = temp5b[2];

                tri[tri_count][0] = temp4b[1];
                tri[tri_count][1] = temp4b[2];
                tri[tri_count][2] = temp3b[1];
                tri_count++;
            }
            else if(far_clip[0] == 2 && far_clip[1] == 0 && far_clip[2] == 1) { // dont split
                tri[tri_count][0] = temp3b[0];
                tri[tri_count][1] = temp3b[1];
                tri[tri_count][2] = temp5b[2];

            }
            else if(far_clip[0] == 3 && far_clip[1] == 2 && far_clip[2] == 1) { // split
                tri[i][0] = temp3b[0];
                tri[i][1] = temp3b[1];
                tri[i][2] = temp4b[2];

                tri[tri_count][0] = temp4b[2];
                tri[tri_count][1] = temp5b[2];
                tri[tri_count][2] = temp3b[0];
                tri_count++;
            }
            else {
                tri_count++;
            }
            
        }
    }
    return tri_count;
}

int clip_line(Vertex v[2]) {
    Vertex center = center_of_polygon(v, 2);
    double x = center.x / center.z;
    double y = center.y / center.z;

    int clip_near = clip_vertex_near(&v[0], &v[1]);
    int clip_far = clip_vertex_far(&v[0], &v[1]);

    if(clip_near == 0 || clip_far == 0) {
        return 0;
    }

    return 1;
}

Vertex calculateForwardVector(double sinYaw, double cosYaw, double sinPitch, double cosPitch) {
    Vertex forward;
    forward.x = cosPitch * sinYaw;
    forward.y = sinPitch;
    forward.z = cosPitch * cosYaw;

    return forward;
}



void scaleLineDepth(Vertex p[2]) {
    p[0].z = vz_to_pz(pz_to_vz(p[0].z)+.02) * .999;
    p[1].z = vz_to_pz(pz_to_vz(p[1].z)+.02) * .999;
}


void drawRubiksCubeTris(Pixel **pixels, RubiksCube *rubiksCube, Vertex rubiksCubePosition, double rubiksCubeSize, Vertex playerPosition, Vertex lightSourcePosition, Matrix rotX, Matrix rotY, Matrix invRotX, Matrix invRotY, bool turn_cube, int axis, int layer) {
    
    for(int i = 0; i < RUBIKS_CUBE_DIM; i++) {
        for(int j = 0; j < RUBIKS_CUBE_DIM; j++) {
            for(int k = 0; k < RUBIKS_CUBE_DIM; k++) {

                Cube cube = rubiksCube->cubies[i][j][k]; 

                for(int f = 0; f < 12; f++) {

                    rubiksCube->cubies[i][j][k].faces[f].isCounterClockwise = false;
                    
                    Vertex v[3] = {cube.faces[f].tri[0], cube.faces[f].tri[1], cube.faces[f].tri[2]};
                    Vertex normal = findNormal(v);
                    int color = cube.faces[f].color;
                    
                    Vertex surfacePosition = normal;
                    surfacePosition.x *= (double)RUBIKS_CUBE_DIM/2.0;
                    surfacePosition.y *= (double)RUBIKS_CUBE_DIM/2.0;
                    surfacePosition.z *= (double)RUBIKS_CUBE_DIM/2.0;

                    bool is_correct_layer = false;
                    if(turn_cube == true) {
                        switch(axis) {
                            case 0:
                                is_correct_layer = (i == layer);
                                break;
                            case 1:
                                is_correct_layer = (j == layer);
                                break;
                            case 2:
                                is_correct_layer = (k == layer);
                                break;
                            default:
                                is_correct_layer = false;
                                break;
                        }
                    }

                    bool is_face_in_turning_layer = isFaceInTurningLayer(v, axis, layer, rubiksCubePosition, rubiksCubeSize);

                    if(isOuterFace(v, rubiksCubeSize, turn_cube, is_correct_layer, is_face_in_turning_layer, axis, layer, *rubiksCube, &color) == false) {
                        continue;
                    }   
        
                    for(int n = 0; n < 3; n++) 
                    {
                        // displace in relation to player
                        v[n].x -= playerPosition.x;
                        v[n].y -= playerPosition.y;
                        v[n].z -= playerPosition.z;

                        // displace for rubiks cube position
                        v[n].x += rubiksCubePosition.x;
                        v[n].y += rubiksCubePosition.y;
                        v[n].z += rubiksCubePosition.z;

                        // camera yaw transform
                        v[n] = vertex_rotate(v[n], rotY);
                        // camera pitch transform
                        v[n] = vertex_rotate(v[n], rotX);
                    }


                    // tri array, usually only 1, at most 2 (if clipped into 2)
                    Vertex tri[2][3];
                    tri[0][0] = v[0];
                    tri[0][1] = v[1];
                    tri[0][2] = v[2];
                    tri[1][0] = v[0];
                    tri[1][1] = v[1];
                    tri[1][2] = v[2];

                    // get clipped triangles (if necessary)
                    int tri_count = clip_tri(v, tri); // will be 0 if all outside of the near plane or far plane

                    // draw triangle(s)
                    for(int t = 0; t < tri_count; t++) {
                        Vertex p[3];
                        p[0] = vertex_project(tri[t][0], projection);
                        p[1] = vertex_project(tri[t][1], projection);
                        p[2] = vertex_project(tri[t][2], projection);

                        Vertex p_view[3];
                        for(int n = 0; n < 3; n++) {
                            p_view[n].x = (((p[n].x + 1.0) * (double)WINDOW_WIDTH) / 2.0);
                            p_view[n].y = (((-p[n].y + 1.0) * (double)WINDOW_HEIGHT) / 2.0);
                        }

                        if(isCounterClockwiseTri(p_view) == 1) {
                            rubiksCube->cubies[i][j][k].faces[f].isCounterClockwise = true;
                            Vertex *oppositeTriAddress;
                            if(f % 2 == 0) {
                                oppositeTriAddress = (Vertex *)rubiksCube->cubies[i][j][k].faces[f+1].tri;
                            }
                            else {
                                oppositeTriAddress = (Vertex *)rubiksCube->cubies[i][j][k].faces[f-1].tri;
                            }
                            // add tri to `pixels`
                            getTriPixels(pixels, p_view, p, normal, surfacePosition, color, lightSourcePosition, playerPosition, invRotX, invRotY, (uint32_t *)(rubiksCube));//(Vertex *)rubiksCube->cubies[i][j][k].faces[f].tri, rubiksCubePosition);
                        }
                    }
                }
            }
        }
    }
}

void drawRubiksCubeLines(Pixel **pixels, RubiksCube *rubiksCube, Vertex rubiksCubePosition, double rubiksCubeSize, Vertex playerPosition, Vertex lightSourcePosition, Matrix rotX, Matrix rotY, Matrix invRotX, Matrix invRotY, bool turn_cube, int axis, int layer) {

    for(int i = 0; i < RUBIKS_CUBE_DIM; i++) {
        for(int j = 0; j < RUBIKS_CUBE_DIM; j++) {
            for(int k = 0; k < RUBIKS_CUBE_DIM; k++) {

                Cube cube = rubiksCube->cubies[i][j][k];

                for(int f = 0; f < 12; f++) {
                    
                    Vertex edges[2][2] = {
                        {cube.faces[f].tri[0], cube.faces[f].tri[1]},
                        {cube.faces[f].tri[1], cube.faces[f].tri[2]},
                    };

                    if(cube.faces[f].isCounterClockwise == false) {
                        continue;
                    }

                    Vertex surfaceNormal = findNormal(cube.faces[f].tri);

                    for(int e = 0; e < 2; e++) {

                        Vertex v[2] = {edges[e][0], edges[e][1]};

                        bool is_correct_layer = false;
                        if(turn_cube == true) {
                            switch(axis) {
                                case 0:
                                    is_correct_layer = (i == layer);
                                    break;
                                case 1:
                                    is_correct_layer = (j == layer);
                                    break;
                                case 2:
                                    is_correct_layer = (k == layer);
                                    break;
                                default:
                                    is_correct_layer = false;
                                    break;
                            }
                        }

                        if(isOuterLine(v, rubiksCubeSize, turn_cube, is_correct_layer, axis, layer, *rubiksCube) == false) {
                            continue;
                        } 

                        for(int n = 0; n < 2; n++) 
                        {
                            // displace in relation to player
                            v[n].x -= playerPosition.x;
                            v[n].y -= playerPosition.y;
                            v[n].z -= playerPosition.z;

                            // displace for rubiks cube position
                            v[n].x += rubiksCubePosition.x;
                            v[n].y += rubiksCubePosition.y;
                            v[n].z += rubiksCubePosition.z;

                            // camera yaw transform
                            v[n] = vertex_rotate(v[n], rotY);
                            // camera pitch transform
                            v[n] = vertex_rotate(v[n], rotX);
                        }

                        // clipping
                        int clip = clip_line(v);

                        if(clip == 1) {
                            // perspective projection transform
                            Vertex p[2];
                            p[0] = vertex_project(v[0], projection);
                            p[1] = vertex_project(v[1], projection);
                            scaleLineDepth(p); // scale depth

                            // view transform
                            Vertex p_view[2];
                            for(int n = 0; n < 2; n++) {
                                p_view[n].x = (((p[n].x + 1.0) * (double)WINDOW_WIDTH) / 2.0);
                                p_view[n].y = (((-p[n].y + 1.0) * (double)WINDOW_HEIGHT) / 2.0);
                            }

                            // add line to `pixels`
                            int color = BLACK;
                            bool isQuadDiag = false;
                            if(e >= 2) {
                                color = cube.faces[f].color;
                                isQuadDiag = true;
                            }
                            getLinePixels(pixels, p_view, p, color, playerPosition, invRotX, invRotY, (uint32_t *)(rubiksCube));
                        }
                    }
                }
            }
        }
    }
}

void drawFloorGridTris(Pixel **pixels, Square *floorGridSquare, Vertex playerPosition, Vertex lightSourcePosition, Matrix rotX, Matrix rotY, Matrix invRotX, Matrix invRotY) {

    for(int f = 0; f < 2; f++) {
            
        Vertex v[3] = {floorGridSquare->faces[f].tri[0], floorGridSquare->faces[f].tri[1], floorGridSquare->faces[f].tri[2]};
        Vertex avg = {
            (v[0].x + v[1].x + v[2].x) / 3.0,
            (v[0].y + v[1].y + v[2].y) / 3.0,
            (v[0].z + v[1].z + v[2].z) / 3.0
        };
        Vertex normal = findNormal(v);
        int color = floorGridSquare->faces[f].color;

        Vertex surfacePosition = {avg.x, avg.y, avg.z};

        for(int n = 0; n < 3; n++) 
        {
            // displace in relation to player
            v[n].x -= playerPosition.x;
            v[n].y -= playerPosition.y;
            v[n].z -= playerPosition.z;

            // camera yaw transform
            v[n] = vertex_rotate(v[n], rotY);
            // camera pitch transform
            v[n] = vertex_rotate(v[n], rotX);
        }

        // tri array
        Vertex tri[4][3];
        for(int t = 0; t < 4; t++) {
            for(int u = 0; u < 3; u++) {
                tri[t][u] = v[u];
            }
        }

        // get clipped triangles (if necessary)
        int tri_count = clip_tri(v, tri); // `tri_count` will be 0 if all outside of the near plane or far plane

        // draw triangle(s)
        for(int t = 0; t < tri_count; t++) {
            Vertex p[3];

            // projection tranform
            p[0] = vertex_project(tri[t][0], projection);
            p[1] = vertex_project(tri[t][1], projection);
            p[2] = vertex_project(tri[t][2], projection);

            Vertex p_view[3];
            for(int n = 0; n < 3; n++) {
                p_view[n].x = (((p[n].x + 1.0) * (double)WINDOW_WIDTH) / 2.0);
                p_view[n].y = (((-p[n].y + 1.0) * (double)WINDOW_HEIGHT) / 2.0);
            }

            if(isCounterClockwiseTri(p_view) == 1) {
                Vertex *oppositeTriAddress;
                if(f == 0) {
                    oppositeTriAddress = (Vertex *)floorGridSquare->faces[f+1].tri;
                }
                else {
                    oppositeTriAddress = (Vertex *)floorGridSquare->faces[f-1].tri;
                }
                // add tri to `pixels`
                getTriPixels(pixels, p_view, p, normal, surfacePosition, color, lightSourcePosition, playerPosition, invRotX, invRotY, (uint32_t *)(floorGridSquare));
            }
        }
    }
}

void drawFloorGridLines(Pixel **pixels, FloorGrid *floorGrid, Square *floorGridSquare, Vertex playerPosition, Vertex lightSourcePosition, Matrix rotX, Matrix rotY, Matrix invRotX, Matrix invRotY) {

    for(int i = 0; i < FLOOR_GRID_DIM; i++) {
        for(int j = 0; j < FLOOR_GRID_DIM; j++) {

            Square square = floorGrid->squares[i][j];

            for(int e = 0; e < 4; e++) {

                Vertex v[2] = {square.edges[e][0], square.edges[e][1]};

                for(int n = 0; n < 2; n++) 
                {
                    // displace in relation to player
                    v[n].x -= playerPosition.x;
                    v[n].y -= playerPosition.y;
                    v[n].z -= playerPosition.z;

                    // camera yaw transform
                    v[n] = vertex_rotate(v[n], rotY);
                    // camera pitch transform
                    v[n] = vertex_rotate(v[n], rotX);
                }

                // clipping
                int clip = clip_line(v);

                if(clip == 1) {
                    // perspective projection transform
                    Vertex p[2];
                    p[0] = vertex_project(v[0], projection);
                    p[1] = vertex_project(v[1], projection);
                    scaleLineDepth(p); // scale depth   

                    // view tranform
                    Vertex p_view[2];
                    for(int n = 0; n < 2; n++) {
                        p_view[n].x = (((p[n].x + 1.0) * (double)WINDOW_WIDTH) / 2.0);
                        p_view[n].y = (((-p[n].y + 1.0) * (double)WINDOW_HEIGHT) / 2.0);
                    }

                    // add line to `pixels`
                    getLinePixels(pixels, p_view, p, BLACK, playerPosition, invRotX, invRotY, (uint32_t *)(floorGridSquare));
                }
            }
        }
    }
}

void drawLightSourceCubeTris(Pixel **pixels, double cube_size, Vertex playerPosition, Vertex lightSourcePosition, Matrix rotX, Matrix rotY, Matrix invRotX, Matrix invRotY) {

    Cube cube = createCube((Vertex){0.0, 0.0, 0.0}, cube_size);

    for(int f = 0; f < 12; f++) {
        
        Vertex v[3] = {cube.faces[f].tri[0], cube.faces[f].tri[1], cube.faces[f].tri[2]};
        Vertex normal = findNormal(v);
        int color = WHITE;

        Vertex surfacePosition = normal;

        normal.x *= -1.0;
        normal.y *= -1.0;
        normal.z *= -1.0;

        surfacePosition.x *= cube_size * 0.5;
        surfacePosition.y *= cube_size * 0.5;
        surfacePosition.z *= cube_size * 0.5;
        surfacePosition.x += lightSourcePosition.x;
        surfacePosition.y += lightSourcePosition.y;
        surfacePosition.z += lightSourcePosition.z;

        for(int n = 0; n < 3; n++) 
        {
            // displace in relation to player
            v[n].x -= playerPosition.x;
            v[n].y -= playerPosition.y;
            v[n].z -= playerPosition.z;

            // displace for light source position
            v[n].x += lightSourcePosition.x;
            v[n].y += lightSourcePosition.y;
            v[n].z += lightSourcePosition.z;

            // camera yaw transform
            v[n] = vertex_rotate(v[n], rotY);
            // camera pitch transform
            v[n] = vertex_rotate(v[n], rotX);
        }

        // tri array, usually only 1, at most 2 (if clipped into 2)
        Vertex tri[2][3];
        tri[0][0] = v[0];
        tri[0][1] = v[1];
        tri[0][2] = v[2];
        tri[1][0] = v[0];
        tri[1][1] = v[1];
        tri[1][2] = v[2];

        // get clipped triangles (if necessary)
        int tri_count = clip_tri(v, tri); // will be 0 if all outside of the near plane or far plane

        // draw triangle(s)
        for(int t = 0; t < tri_count; t++) {
            // perspective projection transform
            Vertex p[3];
            p[0] = vertex_project(tri[t][0], projection);
            p[1] = vertex_project(tri[t][1], projection);
            p[2] = vertex_project(tri[t][2], projection);

            // view transform
            Vertex p_view[3];
            for(int n = 0; n < 3; n++) {
                p_view[n].x = (((p[n].x + 1.0) * (double)WINDOW_WIDTH) / 2.0);
                p_view[n].y = (((-p[n].y + 1.0) * (double)WINDOW_HEIGHT) / 2.0);
            }

            if(isCounterClockwiseTri(p_view) == 1) {
                bool isEvenQuadFace = (f % 2 == 0);
                //Vertex *oppositeTriAddress;
                //if(f % 2 == 0) {
                    //oppositeTriAddress = (Vertex *)cube->faces[f+1].tri;
                //}
                //else {
                    //oppositeTriAddress = (Vertex *)cube->faces[f-1].tri;
                //}
                // add tri to `pixels`
                getTriPixels(pixels, p_view, p, normal, surfacePosition, color, lightSourcePosition, playerPosition, invRotX, invRotY, (uint32_t *)NULL);
            }
        }
    }
}

void drawLightSourceCubeLines(Pixel **pixels, double cube_size, Vertex playerPosition, Vertex lightSourcePosition, Matrix rotX, Matrix rotY, Matrix invRotX, Matrix invRotY) {

    Cube cube = createCube((Vertex){0, 0, 0}, cube_size);

    for(int e = 0; e < 12; e++) {

        Vertex v[2] = {cube.edges[e][0], cube.edges[e][1]};

        for(int n = 0; n < 2; n++) 
        {
            // displace in relation to player
            v[n].x -= playerPosition.x;
            v[n].y -= playerPosition.y;
            v[n].z -= playerPosition.z;

            // displace for light source position
            v[n].x += lightSourcePosition.x;
            v[n].y += lightSourcePosition.y;
            v[n].z += lightSourcePosition.z;

            // camera yaw transform
            v[n] = vertex_rotate(v[n], rotY);
            // camera pitch transform
            v[n] = vertex_rotate(v[n], rotX);
        }

        // clipping
        int clip = clip_line(v);

        if(clip == 1) {
            // perspective projection transform
            Vertex p[2];
            p[0] = vertex_project(v[0], projection);
            p[1] = vertex_project(v[1], projection);
            scaleLineDepth(p); // scale depth

            // view transform
            Vertex p_view[2];
            for(int n = 0; n < 2; n++) {
                p_view[n].x = (((p[n].x + 1.0) * (double)WINDOW_WIDTH) / 2.0);
                p_view[n].y = (((-p[n].y + 1.0) * (double)WINDOW_HEIGHT) / 2.0);
            }

            // add line to `pixels`
            getLinePixels(pixels, p_view, p, BLACK, playerPosition, invRotX, invRotY, (uint32_t *)NULL);
        }
    }
}
















#define PYRAMINX_EDGE_LENGTH 1.224744871 // 1.0 / sin(acos(1.0 / sqrt(3.0)));  -  this ensures that the apex's Y-value will equal 1.0


typedef struct {
    Vertex vertices[4];
    Vertex edges[6][2];
    Face faces[4];
} Tetrahedron;

/*
Vertex *calculateTetrahedronVertices(double edgeLength) {

    Vertex *v = malloc(sizeof(Vertex) * 4);

    // Calculate the height of the equilateral triangle that forms the base
    double h = sqrt(3.0) / 2.0 * edgeLength;

    // Calculate the height of the tetrahedron from its base to the apex
    double H = sqrt(2.0 / 3.0) * edgeLength;

    // front right of base
    v[0].x = edgeLength / 2.0;
    v[0].y = -h / 3.0;
    v[0].z = 0.0;

    // front left of base
    v[1].x = -edgeLength / 2.0;
    v[1].y = -h / 3.0;
    v[1].z = 0.0;

    // back middle of base
    v[2].x = 0.0;
    v[2].y = 2.0 * h / 3.0;
    v[2].z = 0.0;

    // apex above the base
    v[3].x = 0.0;
    v[3].y = 0.0;
    v[3].z = H;

    return v;
}*/

Tetrahedron createTetrahedron() {

    Vertex v[4];

    double edgeLength = PYRAMINX_EDGE_LENGTH;

    double y = 4.0 / 9.0;

    // Calculate the height of the equilateral triangle that forms the base
    double h = 1.73205080757 / 2.0 * edgeLength; //sqrt(3.0) / 2.0 * edgeLength;

    // Calculate the height of the tetrahedron from its base to the apex
    double H = 0.81649658092 * edgeLength; //sqrt(2.0 / 3.0) * edgeLength;

    // front right of base
    v[0].x = -edgeLength / 2.0;
    v[0].z = h / 3.0;
    v[0].y = -y;

    // front left of base
    v[1].x = edgeLength / 2.0;
    v[1].z = h / 3.0;
    v[1].y = -y;

    // back middle of base
    v[2].x = 0.0;
    v[2].z = (-2.0 * h) / 3.0;
    v[2].y = -y;

    // apex above the base
    v[3].x = 0.0;
    v[3].z = 0.0;
    v[3].y = H - y;

    Tetrahedron tetrahedron = {

        .vertices = {v[0], v[1], v[2], v[3]},

        .edges = {
            {v[0], v[1]},
            {v[1], v[2]},
            {v[2], v[0]},
            {v[0], v[3]},
            {v[1], v[3]},
            {v[2], v[3]}
        },

        .faces[0] = {
            .tri = {v[0], v[1], v[2]},
            .color = YELLOW
        },
        .faces[1] = {
            .tri = {v[0], v[3], v[1]},
            .color = RED
        },
        .faces[2] = {
            .tri = {v[1], v[3], v[2]},
            .color = GREEN
        },
        .faces[3] = {
            .tri = {v[2], v[3], v[0]},
            .color = BLUE
        }
    };

    return tetrahedron;
}

Tetrahedron updateTetrahedronEdgesAndFaces(Tetrahedron tetrahedron) {

    Tetrahedron newTetrahedron = {

        .vertices = {tetrahedron.vertices[0], tetrahedron.vertices[1], tetrahedron.vertices[2], tetrahedron.vertices[3]},

        .edges = {
            {tetrahedron.vertices[0], tetrahedron.vertices[1]},
            {tetrahedron.vertices[1], tetrahedron.vertices[2]},
            {tetrahedron.vertices[2], tetrahedron.vertices[0]},
            {tetrahedron.vertices[0], tetrahedron.vertices[3]},
            {tetrahedron.vertices[1], tetrahedron.vertices[3]},
            {tetrahedron.vertices[2], tetrahedron.vertices[3]}
        },

        .faces[0] = {
            .tri = {tetrahedron.vertices[0], tetrahedron.vertices[1], tetrahedron.vertices[2]},
            .color = YELLOW
        },
        .faces[1] = {
            .tri = {tetrahedron.vertices[0], tetrahedron.vertices[3], tetrahedron.vertices[1]},
            .color = RED
        },
        .faces[2] = {
            .tri = {tetrahedron.vertices[1], tetrahedron.vertices[3], tetrahedron.vertices[2]},
            .color = GREEN
        },
        .faces[3] = {
            .tri = {tetrahedron.vertices[2], tetrahedron.vertices[3], tetrahedron.vertices[0]},
            .color = BLUE
        }
    };

    return newTetrahedron;
}





typedef struct {
    Vertex vertices[6];
    Vertex edges[12][2];
    Face faces[8];
} Octahedron;




Octahedron createOctahedron() {

    Vertex v[6];

    double edgeLength = PYRAMINX_EDGE_LENGTH;

    double y = 4.0 / 9.0;

    // Calculate the height of the equilateral triangle that forms the base
    double h = (1.73205080757 / 2.0) * edgeLength; //sqrt(3.0) / 2.0 * edgeLength;

    // Calculate the height of the tetrahedron from its base to the apex
    double H = 0.81649658092 * edgeLength; //sqrt(2.0 / 3.0) * edgeLength;

    // front left of base
    v[0].x = -edgeLength / 2.0;
    v[0].z = (-h / 3.0) - (h / 3.0);
    v[0].y = -y;

    // front right of base
    v[1].x = edgeLength / 2.0;
    v[1].z = (-h / 3.0) - (h / 3.0);
    v[1].y = -y;

    // back middle of base
    v[2].x = 0.0;
    v[2].z = ((2.0 * h) / 3.0) - (h / 3.0);
    v[2].y = -y;



    // front left of top
    v[3].x = -edgeLength / 2.0;
    v[3].z = (h / 3.0) - (h / 3.0);
    v[3].y = H - y;

    // front left of top
    v[4].x = edgeLength / 2.0;
    v[4].z = (h / 3.0) - (h / 3.0);
    v[4].y = H - y;

    // back middle of top
    v[5].x = 0.0;
    v[5].z = ((-2.0 * h) / 3.0) - (h / 3.0);
    v[5].y = H - y;

    



    Octahedron octahedron = {
        .vertices = {v[0], v[1], v[2], v[3], v[4], v[5]},
        .edges = {
            {v[0], v[1]},
            {v[1], v[2]},
            {v[2], v[0]},
            {v[3], v[4]},
            {v[4], v[5]},
            {v[5], v[3]},
            {v[0], v[3]}
        },
    };

    return octahedron;
}






Octahedron updateOctahedronEdgesAndFaces(Octahedron octahedron) {

    Octahedron newOctahedron = {

        .vertices = {octahedron.vertices[0], octahedron.vertices[1], octahedron.vertices[2], octahedron.vertices[3], octahedron.vertices[4], octahedron.vertices[5]},

        .edges = {
            {octahedron.vertices[0], octahedron.vertices[1]},
            {octahedron.vertices[1], octahedron.vertices[2]},
            {octahedron.vertices[2], octahedron.vertices[0]},
            {octahedron.vertices[3], octahedron.vertices[4]},
            {octahedron.vertices[4], octahedron.vertices[5]},
            {octahedron.vertices[5], octahedron.vertices[3]}
        },

        .faces[0] = {
            .tri = {octahedron.vertices[0], octahedron.vertices[1], octahedron.vertices[2]},
            .color = YELLOW
        },
        .faces[1] = {
            .tri = {octahedron.vertices[5], octahedron.vertices[4], octahedron.vertices[3]},
            .color = YELLOW
        },
        .faces[2] = {
            .tri = {octahedron.vertices[0], octahedron.vertices[3], octahedron.vertices[4]},
            .color = RED
        },
        .faces[3] = {
            .tri = {octahedron.vertices[2], octahedron.vertices[1], octahedron.vertices[5]},
            .color = RED
        },
        .faces[4] = {
            .tri = {octahedron.vertices[1], octahedron.vertices[4], octahedron.vertices[5]},
            .color = GREEN
        },
        .faces[5] = {
            .tri = {octahedron.vertices[3], octahedron.vertices[0], octahedron.vertices[2]},
            .color = GREEN
        },
        .faces[6] = {
            .tri = {octahedron.vertices[4], octahedron.vertices[1], octahedron.vertices[0]},
            .color = BLUE
        },
        .faces[7] = {
            .tri = {octahedron.vertices[2], octahedron.vertices[5], octahedron.vertices[3]},
            .color = BLUE
        }
    };

    return newOctahedron;
}




#define TOTAL_TETRAHEDRA 10
#define TOTAL_OCTAHEDRA 4

typedef struct {
    //Vertex t0_position;
    Tetrahedron tetrahedra[10];
    Octahedron octahedra[4];
    double layer_angle[4][3];
    Vertex axis_vectors[8];
} Pyraminx;

Pyraminx createPyraminx(Pyraminx pyraminx) {
    double x_offset = 0.0;
    double z_offset = 0.0;
    double y_offset = 0.0;


    Matrix rotX = rotation_x(degrees_to_radians(70.5288));
    Matrix invRotX = rotation_x(degrees_to_radians(-70.5288));
    Matrix rotY = rotation_y(degrees_to_radians(-180.0));

    Matrix rotX120 = rotation_x(degrees_to_radians(120.0));

    for(int i = 0; i < 1; i++) {
        x_offset = 0.0;
        y_offset = 0.0;
        z_offset = 0.0;
        Vertex tetrahedronPosition = {0, 0, 0};

        Tetrahedron tetrahedron = createTetrahedron();

        for(int j = 0; j < 4; j++) {
            tetrahedron.vertices[j].x += x_offset;
            tetrahedron.vertices[j].y += y_offset;
            tetrahedron.vertices[j].z += z_offset;
        }

        Tetrahedron updatedTetrahedron = updateTetrahedronEdgesAndFaces(tetrahedron);
        pyraminx.tetrahedra[i] = updatedTetrahedron;
    }  
    
    for(int i = 1; i < 3; i++) {
        x_offset = 0.0;
        y_offset = 0.0;
        z_offset = 0.0;
        Vertex tetrahedronPosition = {0, 0, 0};

        Tetrahedron tetrahedron;
        tetrahedron.vertices[0].x = pyraminx.tetrahedra[i-1].vertices[1].x;
        tetrahedron.vertices[1].x = tetrahedron.vertices[0].x + PYRAMINX_EDGE_LENGTH;
        tetrahedron.vertices[2].x = tetrahedron.vertices[0].x + (PYRAMINX_EDGE_LENGTH / 2.0);
        tetrahedron.vertices[3].x = tetrahedron.vertices[0].x + (PYRAMINX_EDGE_LENGTH / 2.0);

        for(int j = 0; j < 4; j++) {
            tetrahedron.vertices[j].z = pyraminx.tetrahedra[i-1].vertices[j].z;
            tetrahedron.vertices[j].y = pyraminx.tetrahedra[i-1].vertices[j].y;
        }

        for(int j = 0; j < 4; j++) {
            tetrahedron.vertices[j].x += x_offset;
            tetrahedron.vertices[j].y += y_offset;
            tetrahedron.vertices[j].z += z_offset;
        }

        Tetrahedron updatedTetrahedron = updateTetrahedronEdgesAndFaces(tetrahedron);
        pyraminx.tetrahedra[i] = updatedTetrahedron;
    }  


    for(int i = 3; i < 5; i++) {
        x_offset = PYRAMINX_EDGE_LENGTH / 2.0;
        y_offset = 0.0;
        z_offset = -0.5 * 1.73205080757 * PYRAMINX_EDGE_LENGTH;
        Vertex tetrahedronPosition = {0, 0, 0};

        Tetrahedron tetrahedron;
        tetrahedron.vertices[0].x = pyraminx.tetrahedra[i-3].vertices[2].x;
        tetrahedron.vertices[0].y = pyraminx.tetrahedra[i-3].vertices[2].y;
        tetrahedron.vertices[0].z = pyraminx.tetrahedra[i-3].vertices[2].z;

        tetrahedron.vertices[1].x = pyraminx.tetrahedra[i-2].vertices[2].x;
        tetrahedron.vertices[1].y = pyraminx.tetrahedra[i-2].vertices[2].y;
        tetrahedron.vertices[1].z = pyraminx.tetrahedra[i-2].vertices[2].z;

        tetrahedron.vertices[2].x = (tetrahedron.vertices[0].x + tetrahedron.vertices[1].x) / 2.0;
        tetrahedron.vertices[2].y = tetrahedron.vertices[0].y;
        tetrahedron.vertices[2].z = tetrahedron.vertices[0].z - fabs(tetrahedron.vertices[0].z - pyraminx.tetrahedra[i-3].vertices[0].z);

        tetrahedron.vertices[3].x = (tetrahedron.vertices[0].x + tetrahedron.vertices[1].x + tetrahedron.vertices[2].x) / 3.0;
        tetrahedron.vertices[3].y = pyraminx.tetrahedra[0].vertices[3].y;
        tetrahedron.vertices[3].z = (tetrahedron.vertices[0].z + tetrahedron.vertices[1].z + tetrahedron.vertices[2].z) / 3.0;


        for(int j = 0; j < 4; j++) {
            //tetrahedron.vertices[j].x += x_offset;
            //tetrahedron.vertices[j].y += y_offset;
            //tetrahedron.vertices[j].z += z_offset;
        }

        Tetrahedron updatedTetrahedron = updateTetrahedronEdgesAndFaces(tetrahedron);
        pyraminx.tetrahedra[i] = updatedTetrahedron;
    }  



    for(int i = 5; i < 6; i++) {
        x_offset = PYRAMINX_EDGE_LENGTH / 2.0;
        y_offset = 0.0;
        z_offset = -0.5 * 1.73205080757 * PYRAMINX_EDGE_LENGTH;
        Vertex tetrahedronPosition = {0, 0, 0};

        Tetrahedron tetrahedron;
        tetrahedron.vertices[0].x = pyraminx.tetrahedra[i-2].vertices[2].x;
        tetrahedron.vertices[0].y = pyraminx.tetrahedra[i-2].vertices[2].y;
        tetrahedron.vertices[0].z = pyraminx.tetrahedra[i-2].vertices[2].z;

        tetrahedron.vertices[1].x = pyraminx.tetrahedra[i-1].vertices[2].x;
        tetrahedron.vertices[1].y = pyraminx.tetrahedra[i-1].vertices[2].y;
        tetrahedron.vertices[1].z = pyraminx.tetrahedra[i-1].vertices[2].z;

        tetrahedron.vertices[2].x = (tetrahedron.vertices[0].x + tetrahedron.vertices[1].x) / 2.0;
        tetrahedron.vertices[2].y = tetrahedron.vertices[0].y;
        tetrahedron.vertices[2].z = tetrahedron.vertices[0].z - fabs(tetrahedron.vertices[0].z - pyraminx.tetrahedra[i-2].vertices[0].z);

        tetrahedron.vertices[3].x = (tetrahedron.vertices[0].x + tetrahedron.vertices[1].x + tetrahedron.vertices[2].x) / 3.0;
        tetrahedron.vertices[3].y = pyraminx.tetrahedra[0].vertices[3].y;
        tetrahedron.vertices[3].z = (tetrahedron.vertices[0].z + tetrahedron.vertices[1].z + tetrahedron.vertices[2].z) / 3.0;


        for(int j = 0; j < 4; j++) {
            //tetrahedron.vertices[j].x += x_offset;
            //tetrahedron.vertices[j].y += y_offset;
            //tetrahedron.vertices[j].z += z_offset;
        }

        Tetrahedron updatedTetrahedron = updateTetrahedronEdgesAndFaces(tetrahedron);
        pyraminx.tetrahedra[i] = updatedTetrahedron;
    }  



    for(int i = 6; i < 8; i++) {

        Tetrahedron tetrahedron;
        tetrahedron.vertices[0].x = pyraminx.tetrahedra[i-6].vertices[3].x;
        tetrahedron.vertices[0].y = pyraminx.tetrahedra[i-6].vertices[3].y;
        tetrahedron.vertices[0].z = pyraminx.tetrahedra[i-6].vertices[3].z;

        tetrahedron.vertices[1].x = pyraminx.tetrahedra[i-5].vertices[3].x;
        tetrahedron.vertices[1].y = pyraminx.tetrahedra[i-5].vertices[3].y;
        tetrahedron.vertices[1].z = pyraminx.tetrahedra[i-5].vertices[3].z;

        tetrahedron.vertices[2].x = pyraminx.tetrahedra[i-3].vertices[3].x;
        tetrahedron.vertices[2].y = pyraminx.tetrahedra[i-3].vertices[3].y;
        tetrahedron.vertices[2].z = pyraminx.tetrahedra[i-3].vertices[3].z;

        tetrahedron.vertices[3].x = (tetrahedron.vertices[0].x + tetrahedron.vertices[1].x + tetrahedron.vertices[2].x) / 3.0;
        tetrahedron.vertices[3].y = pyraminx.tetrahedra[0].vertices[3].y + 1.0;
        tetrahedron.vertices[3].z = (tetrahedron.vertices[0].z + tetrahedron.vertices[1].z + tetrahedron.vertices[2].z) / 3.0;

        Tetrahedron updatedTetrahedron = updateTetrahedronEdgesAndFaces(tetrahedron);
        pyraminx.tetrahedra[i] = updatedTetrahedron;
    }  

    for(int i = 8; i < 9; i++) {

        Tetrahedron tetrahedron;
        tetrahedron.vertices[0].x = pyraminx.tetrahedra[3].vertices[3].x;
        tetrahedron.vertices[0].y = pyraminx.tetrahedra[3].vertices[3].y;
        tetrahedron.vertices[0].z = pyraminx.tetrahedra[3].vertices[3].z;

        tetrahedron.vertices[1].x = pyraminx.tetrahedra[4].vertices[3].x;
        tetrahedron.vertices[1].y = pyraminx.tetrahedra[4].vertices[3].y;
        tetrahedron.vertices[1].z = pyraminx.tetrahedra[4].vertices[3].z;

        tetrahedron.vertices[2].x = pyraminx.tetrahedra[5].vertices[3].x;
        tetrahedron.vertices[2].y = pyraminx.tetrahedra[5].vertices[3].y;
        tetrahedron.vertices[2].z = pyraminx.tetrahedra[5].vertices[3].z;

        tetrahedron.vertices[3].x = (tetrahedron.vertices[0].x + tetrahedron.vertices[1].x + tetrahedron.vertices[2].x) / 3.0;
        tetrahedron.vertices[3].y = pyraminx.tetrahedra[0].vertices[3].y + 1.0;
        tetrahedron.vertices[3].z = (tetrahedron.vertices[0].z + tetrahedron.vertices[1].z + tetrahedron.vertices[2].z) / 3.0;

        Tetrahedron updatedTetrahedron = updateTetrahedronEdgesAndFaces(tetrahedron);
        pyraminx.tetrahedra[i] = updatedTetrahedron;
    }  

    for(int i = 9; i < 10; i++) {

        Tetrahedron tetrahedron;
        tetrahedron.vertices[0].x = pyraminx.tetrahedra[6].vertices[3].x;
        tetrahedron.vertices[0].y = pyraminx.tetrahedra[6].vertices[3].y;
        tetrahedron.vertices[0].z = pyraminx.tetrahedra[6].vertices[3].z;

        tetrahedron.vertices[1].x = pyraminx.tetrahedra[7].vertices[3].x;
        tetrahedron.vertices[1].y = pyraminx.tetrahedra[7].vertices[3].y;
        tetrahedron.vertices[1].z = pyraminx.tetrahedra[7].vertices[3].z;

        tetrahedron.vertices[2].x = pyraminx.tetrahedra[8].vertices[3].x;
        tetrahedron.vertices[2].y = pyraminx.tetrahedra[8].vertices[3].y;
        tetrahedron.vertices[2].z = pyraminx.tetrahedra[8].vertices[3].z;

        tetrahedron.vertices[3].x = (tetrahedron.vertices[0].x + tetrahedron.vertices[1].x + tetrahedron.vertices[2].x) / 3.0;
        tetrahedron.vertices[3].y = pyraminx.tetrahedra[6].vertices[3].y + 1.0;
        tetrahedron.vertices[3].z = (tetrahedron.vertices[0].z + tetrahedron.vertices[1].z + tetrahedron.vertices[2].z) / 3.0;

        Tetrahedron updatedTetrahedron = updateTetrahedronEdgesAndFaces(tetrahedron);
        pyraminx.tetrahedra[i] = updatedTetrahedron;
    }  






    Vertex center;
    center.x = (pyraminx.tetrahedra[0].vertices[0].x + pyraminx.tetrahedra[2].vertices[1].x + pyraminx.tetrahedra[5].vertices[2].x + pyraminx.tetrahedra[9].vertices[3].x) / 4.0;
    center.y = (pyraminx.tetrahedra[0].vertices[0].y + pyraminx.tetrahedra[2].vertices[1].y + pyraminx.tetrahedra[5].vertices[2].y + pyraminx.tetrahedra[9].vertices[3].y) / 4.0;
    center.z = (pyraminx.tetrahedra[0].vertices[0].z + pyraminx.tetrahedra[2].vertices[1].z + pyraminx.tetrahedra[5].vertices[2].z + pyraminx.tetrahedra[9].vertices[3].z) / 4.0;

    

    for(int i = 0; i < TOTAL_TETRAHEDRA; i++) {
        for(int j = 0; j < 4; j++) {
            pyraminx.tetrahedra[i].vertices[j].x -= center.x;
            pyraminx.tetrahedra[i].vertices[j].y -= center.y;
            pyraminx.tetrahedra[i].vertices[j].z -= center.z;
        }
        pyraminx.tetrahedra[i] = updateTetrahedronEdgesAndFaces(pyraminx.tetrahedra[i]);
    }
    Matrix rotate180Z = rotation_z(degrees_to_radians(180.0));
    for(int i = 0; i < TOTAL_TETRAHEDRA; i++) {
        for(int j = 0; j < 4; j++) {
            //pyraminx.tetrahedra[i].vertices[j].x *= -1;
            //pyraminx.tetrahedra[i].vertices[j].y *= -1;
            //pyraminx.tetrahedra[i].vertices[j].z *= -1;
        }
    }

    Vertex tips[4] = {pyraminx.tetrahedra[0].vertices[0], pyraminx.tetrahedra[2].vertices[1], pyraminx.tetrahedra[5].vertices[2], pyraminx.tetrahedra[9].vertices[3]};
    pyraminx.axis_vectors[0] = findNormal((Vertex [3]){tips[0], tips[1], tips[2]});
    pyraminx.axis_vectors[1] = findNormal((Vertex [3]){tips[0], tips[3], tips[1]});
    pyraminx.axis_vectors[2] = findNormal((Vertex [3]){tips[0], tips[2], tips[3]});
    pyraminx.axis_vectors[3] = findNormal((Vertex [3]){tips[2], tips[1], tips[3]});
    pyraminx.axis_vectors[4] = findNormal((Vertex [3]){tips[0], tips[2], tips[1]});
    pyraminx.axis_vectors[5] = findNormal((Vertex [3]){tips[0], tips[1], tips[3]});
    pyraminx.axis_vectors[6] = findNormal((Vertex [3]){tips[0], tips[3], tips[2]});
    pyraminx.axis_vectors[7] = findNormal((Vertex [3]){tips[2], tips[3], tips[1]});

    /*
    printf("%f, %f, %f\n", pyraminx.tetrahedra[0].vertices[0].x, pyraminx.tetrahedra[0].vertices[0].y, pyraminx.tetrahedra[0].vertices[0].z);
    printf("%f, %f, %f\n", pyraminx.tetrahedra[2].vertices[1].x, pyraminx.tetrahedra[2].vertices[1].y, pyraminx.tetrahedra[2].vertices[1].z);
    printf("%f, %f, %f\n", pyraminx.tetrahedra[5].vertices[2].x, pyraminx.tetrahedra[5].vertices[2].y, pyraminx.tetrahedra[5].vertices[2].z);
    printf("%f, %f, %f\n", pyraminx.tetrahedra[9].vertices[3].x, pyraminx.tetrahedra[9].vertices[3].y, pyraminx.tetrahedra[9].vertices[3].z);
    */



















    for(int i = 0; i < 2; i++) {

        Octahedron octahedron;
        octahedron.vertices[0].x = pyraminx.tetrahedra[i].vertices[1].x;
        octahedron.vertices[0].y = pyraminx.tetrahedra[i].vertices[1].y;
        octahedron.vertices[0].z = pyraminx.tetrahedra[i].vertices[1].z;

        octahedron.vertices[1].x = pyraminx.tetrahedra[i+1].vertices[2].x;
        octahedron.vertices[1].y = pyraminx.tetrahedra[i+1].vertices[2].y;
        octahedron.vertices[1].z = pyraminx.tetrahedra[i+1].vertices[2].z;

        octahedron.vertices[2].x = pyraminx.tetrahedra[i].vertices[2].x;
        octahedron.vertices[2].y = pyraminx.tetrahedra[i].vertices[2].y;
        octahedron.vertices[2].z = pyraminx.tetrahedra[i].vertices[2].z;

        octahedron.vertices[3].x = pyraminx.tetrahedra[i].vertices[3].x;
        octahedron.vertices[3].y = pyraminx.tetrahedra[i].vertices[3].y;
        octahedron.vertices[3].z = pyraminx.tetrahedra[i].vertices[3].z;

        octahedron.vertices[4].x = pyraminx.tetrahedra[i+1].vertices[3].x;
        octahedron.vertices[4].y = pyraminx.tetrahedra[i+1].vertices[3].y;
        octahedron.vertices[4].z = pyraminx.tetrahedra[i+1].vertices[3].z;

        octahedron.vertices[5].x = pyraminx.tetrahedra[i+3].vertices[3].x;
        octahedron.vertices[5].y = pyraminx.tetrahedra[i+3].vertices[3].y;
        octahedron.vertices[5].z = pyraminx.tetrahedra[i+3].vertices[3].z;

        Octahedron updatedOctahedron = updateOctahedronEdgesAndFaces(octahedron);
        pyraminx.octahedra[i] = updatedOctahedron;
    }  


    for(int i = 2; i < 3; i++) {

        Octahedron octahedron;
        octahedron.vertices[0].x = pyraminx.tetrahedra[3].vertices[1].x;
        octahedron.vertices[0].y = pyraminx.tetrahedra[3].vertices[1].y;
        octahedron.vertices[0].z = pyraminx.tetrahedra[3].vertices[1].z;

        octahedron.vertices[1].x = pyraminx.tetrahedra[4].vertices[2].x;
        octahedron.vertices[1].y = pyraminx.tetrahedra[4].vertices[2].y;
        octahedron.vertices[1].z = pyraminx.tetrahedra[4].vertices[2].z;

        octahedron.vertices[2].x = pyraminx.tetrahedra[3].vertices[2].x;
        octahedron.vertices[2].y = pyraminx.tetrahedra[3].vertices[2].y;
        octahedron.vertices[2].z = pyraminx.tetrahedra[3].vertices[2].z;

        octahedron.vertices[3].x = pyraminx.tetrahedra[3].vertices[3].x;
        octahedron.vertices[3].y = pyraminx.tetrahedra[3].vertices[3].y;
        octahedron.vertices[3].z = pyraminx.tetrahedra[3].vertices[3].z;

        octahedron.vertices[4].x = pyraminx.tetrahedra[4].vertices[3].x;
        octahedron.vertices[4].y = pyraminx.tetrahedra[4].vertices[3].y;
        octahedron.vertices[4].z = pyraminx.tetrahedra[4].vertices[3].z;

        octahedron.vertices[5].x = pyraminx.tetrahedra[5].vertices[3].x;
        octahedron.vertices[5].y = pyraminx.tetrahedra[5].vertices[3].y;
        octahedron.vertices[5].z = pyraminx.tetrahedra[5].vertices[3].z;

        Octahedron updatedOctahedron = updateOctahedronEdgesAndFaces(octahedron);
        pyraminx.octahedra[i] = updatedOctahedron;
    }  

    for(int i = 3; i < 4; i++) {

        Octahedron octahedron;
        octahedron.vertices[0].x = pyraminx.tetrahedra[6].vertices[1].x;
        octahedron.vertices[0].y = pyraminx.tetrahedra[6].vertices[1].y;
        octahedron.vertices[0].z = pyraminx.tetrahedra[6].vertices[1].z;

        octahedron.vertices[1].x = pyraminx.tetrahedra[7].vertices[2].x;
        octahedron.vertices[1].y = pyraminx.tetrahedra[7].vertices[2].y;
        octahedron.vertices[1].z = pyraminx.tetrahedra[7].vertices[2].z;

        octahedron.vertices[2].x = pyraminx.tetrahedra[6].vertices[2].x;
        octahedron.vertices[2].y = pyraminx.tetrahedra[6].vertices[2].y;
        octahedron.vertices[2].z = pyraminx.tetrahedra[6].vertices[2].z;

        octahedron.vertices[3].x = pyraminx.tetrahedra[6].vertices[3].x;
        octahedron.vertices[3].y = pyraminx.tetrahedra[6].vertices[3].y;
        octahedron.vertices[3].z = pyraminx.tetrahedra[6].vertices[3].z;

        octahedron.vertices[4].x = pyraminx.tetrahedra[7].vertices[3].x;
        octahedron.vertices[4].y = pyraminx.tetrahedra[7].vertices[3].y;
        octahedron.vertices[4].z = pyraminx.tetrahedra[7].vertices[3].z;

        octahedron.vertices[5].x = pyraminx.tetrahedra[8].vertices[3].x;
        octahedron.vertices[5].y = pyraminx.tetrahedra[8].vertices[3].y;
        octahedron.vertices[5].z = pyraminx.tetrahedra[8].vertices[3].z;

        Octahedron updatedOctahedron = updateOctahedronEdgesAndFaces(octahedron);
        pyraminx.octahedra[i] = updatedOctahedron;
    }  

    return pyraminx;
}










/*
void drawTetrahedronTris(Pixel **pixels, Tetrahedron tetrahedron, Vertex tetrahedronPosition, double edgeLength, Vertex playerPosition, Vertex lightSourcePosition, Matrix rotX, Matrix rotY, Matrix invRotX, Matrix invRotY) {

    for(int f = 0; f < 4; f++) {
        
        Vertex v[3] = {tetrahedron.faces[f].tri[0], tetrahedron.faces[f].tri[1], tetrahedron.faces[f].tri[2]};
        Vertex normal = findNormal(v);
        int color = tetrahedron.faces[f].color;
        
        Vertex surfacePosition = normal;
        surfacePosition.x *= 0.136083;
        surfacePosition.y *= 0.136083;
        surfacePosition.z *= 0.136083;

        for(int n = 0; n < 3; n++) 
        {
            // displace in relation to player
            v[n].x -= playerPosition.x;
            v[n].y -= playerPosition.y;
            v[n].z -= playerPosition.z;

            // displace for rubiks cube position
            v[n].x += tetrahedronPosition.x;
            v[n].y += tetrahedronPosition.y;
            v[n].z += tetrahedronPosition.z;

            // camera yaw transform
            v[n] = vertex_rotate(v[n], rotY);
            // camera pitch transform
            v[n] = vertex_rotate(v[n], rotX);
        }


        // tri array, usually only 1, at most 2 (if clipped into 2)
        Vertex tri[2][3];
        tri[0][0] = v[0];
        tri[0][1] = v[1];
        tri[0][2] = v[2];
        tri[1][0] = v[0];
        tri[1][1] = v[1];
        tri[1][2] = v[2];

        // get clipped triangles (if necessary)
        int tri_count = clip_tri(v, tri); // will be 0 if all outside of the near plane or far plane

        // draw triangle(s)
        for(int t = 0; t < tri_count; t++) {
            Vertex p[3];
            p[0] = vertex_project(tri[t][0], projection);
            p[1] = vertex_project(tri[t][1], projection);
            p[2] = vertex_project(tri[t][2], projection);

            // view transform
            SDL_Point points[3];
            view_transform_tri(p, points);

            if(isCounterClockwise(points) == 1) {
                // add tri to `pixels`
                getTriPixels(pixels, points, p, normal, surfacePosition, color, lightSourcePosition, playerPosition, invRotX, invRotY, tetrahedron.faces[f].tri);
            }
        }
    }
}

void drawTetrahedronLines(Pixel **pixels, Tetrahedron tetrahedron, Vertex tetrahedronPosition, double edgeLength, Vertex playerPosition, Matrix rotX, Matrix rotY, Matrix invRotX, Matrix invRotY) {

    for(int e = 0; e < 6; e++) {

        Vertex v[2] = {tetrahedron.edges[e][0], tetrahedron.edges[e][1]};

        double s = 1.004;
        for(int n = 0; n < 2; n++) 
        {
            // scale lines slightly for tetrahedron
            v[n].x *= s;
            v[n].y *= s;
            v[n].z *= s;

            // displace in relation to player
            v[n].x -= playerPosition.x;
            v[n].y -= playerPosition.y;
            v[n].z -= playerPosition.z;

            // displace for tetrahedron position
            v[n].x += tetrahedronPosition.x;
            v[n].y += tetrahedronPosition.y;
            v[n].z += tetrahedronPosition.z;

            // camera yaw transform
            v[n] = vertex_rotate(v[n], rotY);
            // camera pitch transform
            v[n] = vertex_rotate(v[n], rotX);
        }

        // clipping
        int clip = clip_line(v);

        if(clip == 1) {
            // perspective projection transform
            Vertex p[2];
            p[0] = vertex_project(v[0], projection);
            p[1] = vertex_project(v[1], projection);
            p[0].z /= 1.0003;
            p[1].z /= 1.0003;

            // view transform
            SDL_Point points[2];
            view_transform_line(p, points);

            // add line to `pixels`
            getLinePixels(pixels, points, p, BLACK, playerPosition, invRotX, invRotY);
        }
    }
}
*/















































































bool isOuterPyraminxFace(Vertex v_rot[3], bool turn_pyraminx, bool in_turning_layer, int axis, int layer, Pyraminx pyraminx, int *color)
{
    Vertex v[3] = {v_rot[0], v_rot[1], v_rot[2]};
    if(in_turning_layer == true) {
        Vertex axis_vector = pyraminx.axis_vectors[axis];
        double turning_layer_angle = pyraminx.layer_angle[axis][layer];
        for(int i = 0; i < 3; i++) {
            v[i] = vertex_rotate(v[i], rotation_axis(axis_vector, -turning_layer_angle));
        }
    }

    Vertex normal_vector = normalize(findNormal(v));
    Vertex axis_vectors[8];
    for(int i = 0; i < 8; i++) {
        axis_vectors[i] = normalize(pyraminx.axis_vectors[i]);
    }

    double scale = 100.0;
    int face_axis = -1;

    for(int i = 0; i < 8; i++) {
        if((int)(normal_vector.x*scale) == (int)(axis_vectors[i].x*scale) && (int)(normal_vector.y*scale) == (int)(axis_vectors[i].y*scale) && (int)(normal_vector.z*scale) == (int)(axis_vectors[i].z*scale)) {
            face_axis = i;
        }
    }
    if(face_axis == 4) {
        face_axis = 0;
    }
    else if(face_axis == 5) {
        face_axis = 1;
    }
    else if(face_axis == 6) {
        face_axis = 2;
    }
    else if(face_axis == 7) {
        face_axis = 3;
    }
    //printf("%d\n", face_axis);

    double outer_y = -.75;
    double e = 0.05;

    double face_vertex_edge_angle = acos(1.0 / sqrt(3.0)); // approx. 54.7356
    double face_edge_face_angle = acos(1.0 / 3.0); // approx. 70.5288
    double vertex_center_vertex_angle = acos(-1.0 / 3.0); // approx. 109.4712


    Vertex v_new[3] = {v[0], v[1], v[2]};
    if(face_axis == 1) {
        for(int i = 0; i < 3; i++) {
            v_new[i] = vertex_rotate(v_new[i], rotation_x(-vertex_center_vertex_angle));
        }
    }
    else if(face_axis == 2) {
        for(int i = 0; i < 3; i++) {
            v_new[i] = vertex_rotate(v_new[i], rotation_y(degrees_to_radians(-120.0)));
            v_new[i] = vertex_rotate(v_new[i], rotation_x(-vertex_center_vertex_angle));
        }
    }
    else if(face_axis == 3) {
        for(int i = 0; i < 3; i++) {
            v_new[i] = vertex_rotate(v_new[i], rotation_y(degrees_to_radians(120.0)));
            v_new[i] = vertex_rotate(v_new[i], rotation_x(-vertex_center_vertex_angle));
        }
    }




    int count = 0;
    for(int i = 0; i < 3; i++) {
        if(v_new[i].y < outer_y + e) {
            count++;
        }
    }
    if(count == 3) {
        return true;
    }


    if(turn_pyraminx == true) {
        if(axis == 1) {
            for(int i = 0; i < 3; i++) {
                v[i] = vertex_rotate(v[i], rotation_x(-vertex_center_vertex_angle));
            }
        }
        else if(axis == 2) {
            for(int i = 0; i < 3; i++) {
                v[i] = vertex_rotate(v[i], rotation_y(degrees_to_radians(-120.0)));
                v[i] = vertex_rotate(v[i], rotation_x(-vertex_center_vertex_angle));
            }
        }
        else if(axis == 3) {
            for(int i = 0; i < 3; i++) {
                v[i] = vertex_rotate(v[i], rotation_y(degrees_to_radians(120.0)));
                v[i] = vertex_rotate(v[i], rotation_x(-vertex_center_vertex_angle));
            }
        }


        int inner_count = 0;
        double layer01 = -.75 + 1.0;
        double layer12 = -.75 + 2.0;

        if(layer == 0) {
            for(int i = 0; i < 3; i++) {
                if(v[i].y > layer01 - e && v[i].y < layer01 + e) {
                    inner_count++;
                }
            }
        }
        else if(layer == 1) {
            if((v[0].y > layer01 - e && v[0].y < layer01 + e && v[1].y > layer01 - e && v[1].y < layer01 + e && v[2].y > layer01 - e && v[2].y < layer01 + e) || (v[0].y > layer12 - e && v[0].y < layer12 + e && v[1].y > layer12 - e && v[1].y < layer12 + e && v[2].y > layer12 - e && v[2].y < layer12 + e)) {
                inner_count = 3;
            }

        }
        else if(layer == 2) {
            for(int i = 0; i < 3; i++) {
                if(v[i].y > layer12 - e && v[i].y < layer12 + e) {
                    inner_count++;
                }
            }
        }

        if(inner_count == 3) {
            if(color != NULL) {
                *color = BLACK;
            }
            return true;
        }

    }
    return false;
}

bool isOuterPyraminxFaceWithoutBlackFaces(Vertex v_rot[3], bool turn_pyraminx, bool in_turning_layer, int axis, int layer, Pyraminx pyraminx)
{
    Vertex v[3] = {v_rot[0], v_rot[1], v_rot[2]};
    if(in_turning_layer == true) {
        Vertex axis_vector = pyraminx.axis_vectors[axis];
        double turning_layer_angle = pyraminx.layer_angle[axis][layer];
        for(int i = 0; i < 3; i++) {
            v[i] = vertex_rotate(v[i], rotation_axis(axis_vector, -turning_layer_angle));
        }
    }

    Vertex normal_vector = normalize(findNormal(v));
    Vertex axis_vectors[8];
    for(int i = 0; i < 8; i++) {
        axis_vectors[i] = normalize(pyraminx.axis_vectors[i]);
    }

    double scale = 100.0;
    int face_axis = -1;

    for(int i = 0; i < 8; i++) {
        if((int)(normal_vector.x*scale) == (int)(axis_vectors[i].x*scale) && (int)(normal_vector.y*scale) == (int)(axis_vectors[i].y*scale) && (int)(normal_vector.z*scale) == (int)(axis_vectors[i].z*scale)) {
            face_axis = i;
        }
    }
    if(face_axis == 4) {
        face_axis = 0;
    }
    else if(face_axis == 5) {
        face_axis = 1;
    }
    else if(face_axis == 6) {
        face_axis = 2;
    }
    else if(face_axis == 7) {
        face_axis = 3;
    }
    //printf("%d\n", face_axis);

    double outer_y = -.75;
    double e = 0.05;

    double face_vertex_edge_angle = acos(1.0 / sqrt(3.0)); // approx. 54.7356
    double face_edge_face_angle = acos(1.0 / 3.0); // approx. 70.5288
    double vertex_center_vertex_angle = acos(-1.0 / 3.0); // approx. 109.4712


    Vertex v_new[3] = {v[0], v[1], v[2]};
    if(face_axis == 1) {
        for(int i = 0; i < 3; i++) {
            v_new[i] = vertex_rotate(v_new[i], rotation_x(-vertex_center_vertex_angle));
        }
    }
    else if(face_axis == 2) {
        for(int i = 0; i < 3; i++) {
            v_new[i] = vertex_rotate(v_new[i], rotation_y(degrees_to_radians(-120.0)));
            v_new[i] = vertex_rotate(v_new[i], rotation_x(-vertex_center_vertex_angle));
        }
    }
    else if(face_axis == 3) {
        for(int i = 0; i < 3; i++) {
            v_new[i] = vertex_rotate(v_new[i], rotation_y(degrees_to_radians(120.0)));
            v_new[i] = vertex_rotate(v_new[i], rotation_x(-vertex_center_vertex_angle));
        }
    }




    int count = 0;
    for(int i = 0; i < 3; i++) {
        if(v_new[i].y < outer_y + e) {
            count++;
        }
    }
    if(count == 3) {
        return true;
    }

    return false;
}

bool isOrthogonalToPlane(Vertex line[2], Vertex normal) {
    // Create a vector representing the line
    Vertex lineVec;
    lineVec.x = line[1].x - line[0].x;
    lineVec.y = line[1].y - line[0].y;
    lineVec.z = line[1].z - line[0].z;

    // Calculate the dot product with the normal
    double dot = dotProduct(lineVec, normal);

    // Check if the dot product is close to zero
    return (dot < 0.001f && dot > -0.001f);
}


bool isOuterPyraminxLine(Vertex v_rot[2], bool turn_pyraminx, bool in_turning_layer, int axis, int layer, Pyraminx pyraminx)
{
    Vertex v[2] = {v_rot[0], v_rot[1]};
    if(in_turning_layer == true) {
        Vertex axis_vector = pyraminx.axis_vectors[axis];
        double turning_layer_angle = pyraminx.layer_angle[axis][layer];
        for(int i = 0; i < 2; i++) {
            v[i] = vertex_rotate(v[i], rotation_axis(axis_vector, -turning_layer_angle));
        }
    }



    int face_axis = -1;

    Vertex axis_vectors[4];
    for(int i = 0; i < 4; i++) {
        axis_vectors[i] = normalize(pyraminx.axis_vectors[i]);
    }

    for(int i = 0; i < 4; i++) {
        if(isOrthogonalToPlane(v, axis_vectors[i]) == true) {
            face_axis = i;
        }
    }
    if(face_axis == 3) {
        if(v[0].x < 0 && v[1].x < 0) {
            face_axis = 2;
        }
    }
    if(face_axis == 2 || face_axis == 3) {
        if(v[0].z > 0 && v[1].z > 0 && isOrthogonalToPlane(v, axis_vectors[1]) == true) {
            face_axis = 1;
        }
    }
    if(face_axis != 0) {
        if(v[0].y < -.74 && v[1].y < -.74) {
            face_axis = 0;
        }
    }
    //printf("%d\n", face_axis);

    double outer_y = -.75;
    double e = 0.05;

    double face_vertex_edge_angle = acos(1.0 / sqrt(3.0)); // approx. 54.7356
    double face_edge_face_angle = acos(1.0 / 3.0); // approx. 70.5288
    double vertex_center_vertex_angle = acos(-1.0 / 3.0); // approx. 109.4712


    Vertex v_new[2] = {v[0], v[1]};
    if(face_axis == 1) {
        for(int i = 0; i < 2; i++) {
            v_new[i] = vertex_rotate(v_new[i], rotation_x(-vertex_center_vertex_angle));
        }
    }
    else if(face_axis == 2) {
        for(int i = 0; i < 2; i++) {
            v_new[i] = vertex_rotate(v_new[i], rotation_y(degrees_to_radians(-120.0)));
            v_new[i] = vertex_rotate(v_new[i], rotation_x(-vertex_center_vertex_angle));
        }
    }
    else if(face_axis == 3) {
        for(int i = 0; i < 2; i++) {
            v_new[i] = vertex_rotate(v_new[i], rotation_y(degrees_to_radians(120.0)));
            v_new[i] = vertex_rotate(v_new[i], rotation_x(-vertex_center_vertex_angle));
        }
    }




    int count = 0;
    for(int i = 0; i < 2; i++) {
        if(v_new[i].y < outer_y + e) {
            count++;
        }
    }
    if(count == 2) {
        return true;
    }

    return false;
}



void drawPyraminxTris(Pixel **pixels, Pyraminx *pyraminx, Vertex pyraminxPosition, Vertex playerPosition, Vertex lightSourcePosition, Matrix rotX, Matrix rotY, Matrix invRotX, Matrix invRotY, bool turn_pyraminx, int axis, int layer) {

    for(int i = 0; i < TOTAL_TETRAHEDRA; i++) {

        Tetrahedron tetrahedron = pyraminx->tetrahedra[i];
        
        for(int f = 0; f < 4; f++) {

            pyraminx->tetrahedra[i].faces[f].isCounterClockwise = false;
            
            Vertex v[3] = {tetrahedron.faces[f].tri[0], tetrahedron.faces[f].tri[1], tetrahedron.faces[f].tri[2]};
            Vertex normal = findNormal(v);
            int color = tetrahedron.faces[f].color;


            
            Vertex surfacePosition = normal;
            surfacePosition.x *= 0.136083;
            surfacePosition.y *= 0.136083;
            surfacePosition.z *= 0.136083;



            bool is_correct_layer = false;
            if(turn_pyraminx == true) 
            {
                if(axis == 0) 
                {
                    if(layer == 0) {
                        if(i == 0 || i == 1 || i == 2 || i == 3 || i == 4 || i == 5) {
                            is_correct_layer = true;
                        }
                    }
                    else if(layer == 1) {
                        if(i == 6 || i == 7 || i == 8) {
                            is_correct_layer = true;
                        }
                    }
                    else if(layer == 2) {
                        if(i == 9) {
                            is_correct_layer = true;
                        }
                    }
                }

                else if(axis == 1) 
                {
                    if(layer == 0) {
                        if(i == 0 || i == 1 || i == 2 || i == 6 || i == 7 || i == 9) {
                            is_correct_layer = true;
                        }
                    }
                    else if(layer == 1) {
                        if(i == 3 || i == 4 || i == 8) {
                            is_correct_layer = true;
                        }
                    }
                    else if(layer == 2) {
                        if(i == 5) {
                            is_correct_layer = true;
                        }
                    }
                }
                else if(axis == 2) 
                {
                    if(layer == 0) {
                        if(i == 0 || i == 3 || i == 5 || i == 6 || i == 8 || i == 9) {
                            is_correct_layer = true;
                        }
                    }
                    else if(layer == 1) {
                        if(i == 1 || i == 4 || i == 7) {
                            is_correct_layer = true;
                        }
                    }
                    else if(layer == 2) {
                        if(i == 2) {
                            is_correct_layer = true;
                        }
                    }
                }

                else if(axis == 3) 
                {
                    if(layer == 0) {
                        if(i == 2 || i == 4 || i == 5 || i == 7 || i == 8 || i == 9) {
                            is_correct_layer = true;
                        }
                    }
                    else if(layer == 1) {
                        if(i == 1 || i == 3 || i == 6) {
                            is_correct_layer = true;
                        }
                    }
                    else if(layer == 2) {
                        if(i == 0) {
                            is_correct_layer = true;
                        }
                    }
                }
            }

            if(isOuterPyraminxFace(v, turn_pyraminx, is_correct_layer, axis, layer, *pyraminx, &color) == false) {
                continue;
            } 




            for(int n = 0; n < 3; n++) 
            {
                // displace in relation to player
                v[n].x -= playerPosition.x;
                v[n].y -= playerPosition.y;
                v[n].z -= playerPosition.z;

                // displace for object position
                v[n].x += pyraminxPosition.x;
                v[n].y += pyraminxPosition.y;
                v[n].z += pyraminxPosition.z;

                // camera yaw transform
                v[n] = vertex_rotate(v[n], rotY);
                // camera pitch transform
                v[n] = vertex_rotate(v[n], rotX);
            }


            // tri array, usually only 1, at most 2 (if clipped into 2)
            Vertex tri[2][3];
            tri[0][0] = v[0];
            tri[0][1] = v[1];
            tri[0][2] = v[2];
            tri[1][0] = v[0];
            tri[1][1] = v[1];
            tri[1][2] = v[2];

            // get clipped triangles (if necessary)
            int tri_count = clip_tri(v, tri); // will be 0 if all outside of the near plane or far plane

            // draw triangle(s)
            for(int t = 0; t < tri_count; t++) {
                Vertex p[3];
                p[0] = vertex_project(tri[t][0], projection);
                p[1] = vertex_project(tri[t][1], projection);
                p[2] = vertex_project(tri[t][2], projection);

                // view transform
                Vertex p_view[3];
                for(int n = 0; n < 3; n++) {
                    p_view[n].x = (((p[n].x + 1.0) * (double)WINDOW_WIDTH) / 2.0);
                    p_view[n].y = (((-p[n].y + 1.0) * (double)WINDOW_HEIGHT) / 2.0);
                }

                if(isCounterClockwiseTri(p_view) == 1) {
                    pyraminx->tetrahedra[i].faces[f].isCounterClockwise = true;
                    // add tri to `pixels`
                    getTriPixels(pixels, p_view, p, normal, surfacePosition, color, lightSourcePosition, playerPosition, invRotX, invRotY, (uint32_t *)(pyraminx));
                }
            }
        }
    }




    
    for(int i = 0; i < TOTAL_OCTAHEDRA; i++) {

        Octahedron octahedron = pyraminx->octahedra[i];
        
        for(int f = 0; f < 8; f++) {

            pyraminx->octahedra[i].faces[f].isCounterClockwise = false;
            
            Vertex v[3] = {octahedron.faces[f].tri[0], octahedron.faces[f].tri[1], octahedron.faces[f].tri[2]};
            Vertex normal = findNormal(v);
            int color = octahedron.faces[f].color;
            
            Vertex surfacePosition = normal;
            surfacePosition.x *= 0.136083;
            surfacePosition.y *= 0.136083;
            surfacePosition.z *= 0.136083;

            bool is_correct_layer = false;
            if(turn_pyraminx == true) 
            {
                if(axis == 0) 
                {
                    if(layer == 0) {
                        if(i == 0 || i == 1 || i == 2) {
                            is_correct_layer = true;
                        }
                    }
                    else if(layer == 1) {
                        if(i == 3) {
                            is_correct_layer = true;
                        }
                    }
                }

                else if(axis == 1) 
                {
                    if(layer == 0) {
                        if(i == 0 || i == 1 || i == 3) {
                            is_correct_layer = true;
                        }
                    }
                    else if(layer == 1) {
                        if(i == 2) {
                            is_correct_layer = true;
                        }
                    }
                }

                else if(axis == 2) 
                {
                    if(layer == 0) {
                        if(i == 0 || i == 2 || i == 3) {
                            is_correct_layer = true;
                        }
                    }
                    else if(layer == 1) {
                        if(i == 1) {
                            is_correct_layer = true;
                        }
                    }
                }

                else if(axis == 3) 
                {
                    if(layer == 0) {
                        if(i == 1 || i == 2 || i == 3) {
                            is_correct_layer = true;
                        }
                    }
                    else if(layer == 1) {
                        if(i == 0) {
                            is_correct_layer = true;
                        }
                    }
                }
            }

            if(isOuterPyraminxFace(v, turn_pyraminx, is_correct_layer, axis, layer, *pyraminx, &color) == false) {
                continue;
            } 

            for(int n = 0; n < 3; n++) 
            {
                // displace in relation to player
                v[n].x -= playerPosition.x;
                v[n].y -= playerPosition.y;
                v[n].z -= playerPosition.z;

                // displace for object position
                v[n].x += pyraminxPosition.x;
                v[n].y += pyraminxPosition.y;
                v[n].z += pyraminxPosition.z;

                // camera yaw transform
                v[n] = vertex_rotate(v[n], rotY);
                // camera pitch transform
                v[n] = vertex_rotate(v[n], rotX);
            }


            // tri array, usually only 1, at most 2 (if clipped into 2)
            Vertex tri[2][3];
            tri[0][0] = v[0];
            tri[0][1] = v[1];
            tri[0][2] = v[2];
            tri[1][0] = v[0];
            tri[1][1] = v[1];
            tri[1][2] = v[2];

            // get clipped triangles (if necessary)
            int tri_count = clip_tri(v, tri); // will be 0 if all outside of the near plane or far plane

            // draw triangle(s)
            for(int t = 0; t < tri_count; t++) {
                Vertex p[3];
                p[0] = vertex_project(tri[t][0], projection);
                p[1] = vertex_project(tri[t][1], projection);
                p[2] = vertex_project(tri[t][2], projection);

                // view transform
                Vertex p_view[3];
                for(int n = 0; n < 3; n++) {
                    p_view[n].x = (((p[n].x + 1.0) * (double)WINDOW_WIDTH) / 2.0);
                    p_view[n].y = (((-p[n].y + 1.0) * (double)WINDOW_HEIGHT) / 2.0);
                }

                if(isCounterClockwiseTri(p_view) == 1) {
                    pyraminx->octahedra[i].faces[f].isCounterClockwise = true;
                    // add tri to `pixels`
                    getTriPixels(pixels, p_view, p, normal, surfacePosition, color, lightSourcePosition, playerPosition, invRotX, invRotY, (uint32_t *)(pyraminx));
                }
            }
        }
    }



}

void drawPyraminxLines(Pixel **pixels, Pyraminx *pyraminx, Vertex pyraminxPosition, Vertex playerPosition, Vertex lightSourcePosition, Matrix rotX, Matrix rotY, Matrix invRotX, Matrix invRotY, bool turn_pyraminx, int axis, int layer) {

    for(int i = 0; i < TOTAL_TETRAHEDRA; i++) {

        Tetrahedron tetrahedron = pyraminx->tetrahedra[i];

        for(int f = 0; f < 4; f++) {

            if(pyraminx->tetrahedra[i].faces[f].isCounterClockwise == false) {
                continue;
            }

            Vertex edges[3][2] = {
                {tetrahedron.faces[f].tri[0], tetrahedron.faces[f].tri[1]},
                {tetrahedron.faces[f].tri[1], tetrahedron.faces[f].tri[2]},
                {tetrahedron.faces[f].tri[2], tetrahedron.faces[f].tri[0]}
            };

            for(int e = 0; e < 3; e++) {

                Vertex v[2] = {edges[e][0], edges[e][1]};

                bool is_correct_layer = false;
                if(turn_pyraminx == true) 
                {
                    if(axis == 0) 
                    {
                        if(layer == 0) {
                            if(i == 0 || i == 1 || i == 2 || i == 3 || i == 4 || i == 5) {
                                is_correct_layer = true;
                            }
                        }
                        else if(layer == 1) {
                            if(i == 6 || i == 7 || i == 8) {
                                is_correct_layer = true;
                            }
                        }
                        else if(layer == 2) {
                            if(i == 9) {
                                is_correct_layer = true;
                            }
                        }
                    }

                    else if(axis == 1) 
                    {
                        if(layer == 0) {
                            if(i == 0 || i == 1 || i == 2 || i == 6 || i == 7 || i == 9) {
                                is_correct_layer = true;
                            }
                        }
                        else if(layer == 1) {
                            if(i == 3 || i == 4 || i == 8) {
                                is_correct_layer = true;
                            }
                        }
                        else if(layer == 2) {
                            if(i == 5) {
                                is_correct_layer = true;
                            }
                        }
                    }
                    else if(axis == 2) 
                    {
                        if(layer == 0) {
                            if(i == 0 || i == 3 || i == 5 || i == 6 || i == 8 || i == 9) {
                                is_correct_layer = true;
                            }
                        }
                        else if(layer == 1) {
                            if(i == 1 || i == 4 || i == 7) {
                                is_correct_layer = true;
                            }
                        }
                        else if(layer == 2) {
                            if(i == 2) {
                                is_correct_layer = true;
                            }
                        }
                    }

                    else if(axis == 3) 
                    {
                        if(layer == 0) {
                            if(i == 2 || i == 4 || i == 5 || i == 7 || i == 8 || i == 9) {
                                is_correct_layer = true;
                            }
                        }
                        else if(layer == 1) {
                            if(i == 1 || i == 3 || i == 6) {
                                is_correct_layer = true;
                            }
                        }
                        else if(layer == 2) {
                            if(i == 0) {
                                is_correct_layer = true;
                            }
                        }
                    }
                }

                if(isOuterPyraminxLine(v, turn_pyraminx, is_correct_layer, axis, layer, *pyraminx) == false) {
                    continue;
                } 

                for(int n = 0; n < 2; n++) 
                {

                    // displace in relation to player
                    v[n].x -= playerPosition.x;
                    v[n].y -= playerPosition.y;
                    v[n].z -= playerPosition.z;

                    // displace for object position
                    v[n].x += pyraminxPosition.x;
                    v[n].y += pyraminxPosition.y;
                    v[n].z += pyraminxPosition.z;

                    // camera yaw transform
                    v[n] = vertex_rotate(v[n], rotY);
                    // camera pitch transform
                    v[n] = vertex_rotate(v[n], rotX);
                }

                // clipping
                int clip = clip_line(v);

                if(clip == 1) {
                    // perspective projection transform
                    Vertex p[2];
                    p[0] = vertex_project(v[0], projection);
                    p[1] = vertex_project(v[1], projection);
                    scaleLineDepth(p); // scale depth

                    // view transform
                    Vertex p_view[2];
                    for(int n = 0; n < 2; n++) {
                        p_view[n].x = (((p[n].x + 1.0) * (double)WINDOW_WIDTH) / 2.0);
                        p_view[n].y = (((-p[n].y + 1.0) * (double)WINDOW_HEIGHT) / 2.0);
                    }

                    // add line to `pixels`
                    getLinePixels(pixels, p_view, p, BLACK, playerPosition, invRotX, invRotY, (uint32_t *)(pyraminx));
                }
            }
        }
    }



    for(int i = 0; i < TOTAL_OCTAHEDRA; i++) {

        Octahedron octahedron = pyraminx->octahedra[i];

        for(int f = 0; f < 8; f++) {

            if(pyraminx->octahedra[i].faces[f].isCounterClockwise == false) {
                continue;
            }

            Vertex edges[3][2] = {
                {octahedron.faces[f].tri[0], octahedron.faces[f].tri[1]},
                {octahedron.faces[f].tri[1], octahedron.faces[f].tri[2]},
                {octahedron.faces[f].tri[2], octahedron.faces[f].tri[0]}
            };

            for(int e = 0; e < 3; e++) {

                Vertex v[2] = {edges[e][0], edges[e][1]};

                bool is_correct_layer = false;
                if(turn_pyraminx == true) 
                {
                    if(axis == 0) 
                    {
                        if(layer == 0) {
                            if(i == 0 || i == 1 || i == 2) {
                                is_correct_layer = true;
                            }
                        }
                        else if(layer == 1) {
                            if(i == 3) {
                                is_correct_layer = true;
                            }
                        }
                    }

                    else if(axis == 1) 
                    {
                        if(layer == 0) {
                            if(i == 0 || i == 1 || i == 3) {
                                is_correct_layer = true;
                            }
                        }
                        else if(layer == 1) {
                            if(i == 2) {
                                is_correct_layer = true;
                            }
                        }
                    }

                    else if(axis == 2) 
                    {
                        if(layer == 0) {
                            if(i == 0 || i == 2 || i == 3) {
                                is_correct_layer = true;
                            }
                        }
                        else if(layer == 1) {
                            if(i == 1) {
                                is_correct_layer = true;
                            }
                        }
                    }

                    else if(axis == 3) 
                    {
                        if(layer == 0) {
                            if(i == 1 || i == 2 || i == 3) {
                                is_correct_layer = true;
                            }
                        }
                        else if(layer == 1) {
                            if(i == 0) {
                                is_correct_layer = true;
                            }
                        }
                    }
                }

                if(isOuterPyraminxLine(v, turn_pyraminx, is_correct_layer, axis, layer, *pyraminx) == false) {
                    continue;
                } 

                for(int n = 0; n < 2; n++) 
                {
                    // displace in relation to player
                    v[n].x -= playerPosition.x;
                    v[n].y -= playerPosition.y;
                    v[n].z -= playerPosition.z;

                    // displace for object position
                    v[n].x += pyraminxPosition.x;
                    v[n].y += pyraminxPosition.y;
                    v[n].z += pyraminxPosition.z;

                    // camera yaw transform
                    v[n] = vertex_rotate(v[n], rotY);
                    // camera pitch transform
                    v[n] = vertex_rotate(v[n], rotX);
                }

                // clipping
                int clip = clip_line(v);

                if(clip == 1) {
                    // perspective projection transform
                    Vertex p[2];
                    p[0] = vertex_project(v[0], projection);
                    p[1] = vertex_project(v[1], projection);
                    scaleLineDepth(p); // scale depth

                    // view transform
                    Vertex p_view[2];
                    for(int n = 0; n < 2; n++) {
                        p_view[n].x = (((p[n].x + 1.0) * (double)WINDOW_WIDTH) / 2.0);
                        p_view[n].y = (((-p[n].y + 1.0) * (double)WINDOW_HEIGHT) / 2.0);
                    }

                    // add line to `pixels`
                    getLinePixels(pixels, p_view, p, BLACK, playerPosition, invRotX, invRotY, (uint32_t *)(pyraminx));
                }
            }
        }
    }
}













void updatePyraminxIndexing(Pyraminx *pyraminx, int axis, int layer, int direction) {
    Pyraminx temp = *pyraminx;

    if(layer == 0) {
        if(axis == 0) {
            if(direction == 1) {
                pyraminx->tetrahedra[0] = temp.tetrahedra[5];
                pyraminx->tetrahedra[1] = temp.tetrahedra[3];
                pyraminx->tetrahedra[2] = temp.tetrahedra[0];
                pyraminx->tetrahedra[3] = temp.tetrahedra[4];
                pyraminx->tetrahedra[4] = temp.tetrahedra[1];
                pyraminx->tetrahedra[5] = temp.tetrahedra[2];
                pyraminx->octahedra[0] = temp.octahedra[2];
                pyraminx->octahedra[1] = temp.octahedra[0];
                pyraminx->octahedra[2] = temp.octahedra[1];
            }
            if(direction == -1) {
                pyraminx->tetrahedra[0] = temp.tetrahedra[2];
                pyraminx->tetrahedra[1] = temp.tetrahedra[4];
                pyraminx->tetrahedra[2] = temp.tetrahedra[5];
                pyraminx->tetrahedra[3] = temp.tetrahedra[1];
                pyraminx->tetrahedra[4] = temp.tetrahedra[3];
                pyraminx->tetrahedra[5] = temp.tetrahedra[0];
                pyraminx->octahedra[0] = temp.octahedra[1];
                pyraminx->octahedra[1] = temp.octahedra[2];
                pyraminx->octahedra[2] = temp.octahedra[0];
            }
        }
        if(axis == 1) {
            if(direction == 1) {
                pyraminx->tetrahedra[0] = temp.tetrahedra[2];
                pyraminx->tetrahedra[1] = temp.tetrahedra[7];
                pyraminx->tetrahedra[2] = temp.tetrahedra[9];
                pyraminx->tetrahedra[6] = temp.tetrahedra[1];
                pyraminx->tetrahedra[7] = temp.tetrahedra[6];
                pyraminx->tetrahedra[9] = temp.tetrahedra[0];
                pyraminx->octahedra[0] = temp.octahedra[1];
                pyraminx->octahedra[1] = temp.octahedra[3];
                pyraminx->octahedra[3] = temp.octahedra[0];
            }
            if(direction == -1) {
                pyraminx->tetrahedra[0] = temp.tetrahedra[9];
                pyraminx->tetrahedra[1] = temp.tetrahedra[6];
                pyraminx->tetrahedra[2] = temp.tetrahedra[0];
                pyraminx->tetrahedra[6] = temp.tetrahedra[7];
                pyraminx->tetrahedra[7] = temp.tetrahedra[1];
                pyraminx->tetrahedra[9] = temp.tetrahedra[2];
                pyraminx->octahedra[0] = temp.octahedra[3];
                pyraminx->octahedra[1] = temp.octahedra[0];
                pyraminx->octahedra[3] = temp.octahedra[1];
            }
        }
        if(axis == 2) {
            if(direction == 1) {
                pyraminx->tetrahedra[0] = temp.tetrahedra[9];
                pyraminx->tetrahedra[3] = temp.tetrahedra[6];
                pyraminx->tetrahedra[5] = temp.tetrahedra[0];
                pyraminx->tetrahedra[6] = temp.tetrahedra[8];
                pyraminx->tetrahedra[8] = temp.tetrahedra[3];
                pyraminx->tetrahedra[9] = temp.tetrahedra[5];
                pyraminx->octahedra[0] = temp.octahedra[3];
                pyraminx->octahedra[2] = temp.octahedra[0];
                pyraminx->octahedra[3] = temp.octahedra[2];
            }
            if(direction == -1) {
                pyraminx->tetrahedra[0] = temp.tetrahedra[5];
                pyraminx->tetrahedra[3] = temp.tetrahedra[8];
                pyraminx->tetrahedra[5] = temp.tetrahedra[9];
                pyraminx->tetrahedra[6] = temp.tetrahedra[3];
                pyraminx->tetrahedra[8] = temp.tetrahedra[6];
                pyraminx->tetrahedra[9] = temp.tetrahedra[0];
                pyraminx->octahedra[0] = temp.octahedra[2];
                pyraminx->octahedra[2] = temp.octahedra[3];
                pyraminx->octahedra[3] = temp.octahedra[0];
            }
        }
        if(axis == 3) {
            if(direction == 1) {
                pyraminx->tetrahedra[2] = temp.tetrahedra[5];
                pyraminx->tetrahedra[4] = temp.tetrahedra[8];
                pyraminx->tetrahedra[5] = temp.tetrahedra[9];
                pyraminx->tetrahedra[7] = temp.tetrahedra[4];
                pyraminx->tetrahedra[8] = temp.tetrahedra[7];
                pyraminx->tetrahedra[9] = temp.tetrahedra[2];
                pyraminx->octahedra[1] = temp.octahedra[2];
                pyraminx->octahedra[2] = temp.octahedra[3];
                pyraminx->octahedra[3] = temp.octahedra[1];
            }
            if(direction == -1) {
                pyraminx->tetrahedra[2] = temp.tetrahedra[9];
                pyraminx->tetrahedra[4] = temp.tetrahedra[7];
                pyraminx->tetrahedra[5] = temp.tetrahedra[2];
                pyraminx->tetrahedra[7] = temp.tetrahedra[8];
                pyraminx->tetrahedra[8] = temp.tetrahedra[4];
                pyraminx->tetrahedra[9] = temp.tetrahedra[5];
                pyraminx->octahedra[1] = temp.octahedra[3];
                pyraminx->octahedra[2] = temp.octahedra[1];
                pyraminx->octahedra[3] = temp.octahedra[2];
            }
        }
    }

    if(layer == 1) {
        if(axis == 0) {
            if(direction == 1) {
                pyraminx->tetrahedra[6] = temp.tetrahedra[8];
                pyraminx->tetrahedra[7] = temp.tetrahedra[6];
                pyraminx->tetrahedra[8] = temp.tetrahedra[7];
            }
            if(direction == -1) {
                pyraminx->tetrahedra[6] = temp.tetrahedra[7];
                pyraminx->tetrahedra[7] = temp.tetrahedra[8];
                pyraminx->tetrahedra[8] = temp.tetrahedra[6];
            }
        }
        if(axis == 1) {
            if(direction == 1) {
                pyraminx->tetrahedra[3] = temp.tetrahedra[4];
                pyraminx->tetrahedra[4] = temp.tetrahedra[8];
                pyraminx->tetrahedra[8] = temp.tetrahedra[3];
            }
            if(direction == -1) {
                pyraminx->tetrahedra[3] = temp.tetrahedra[8];
                pyraminx->tetrahedra[4] = temp.tetrahedra[3];
                pyraminx->tetrahedra[8] = temp.tetrahedra[4];
            }
        }
        if(axis == 2) {
            if(direction == 1) {
                pyraminx->tetrahedra[1] = temp.tetrahedra[7];
                pyraminx->tetrahedra[4] = temp.tetrahedra[1];
                pyraminx->tetrahedra[7] = temp.tetrahedra[4];
            }
            if(direction == -1) {
                pyraminx->tetrahedra[1] = temp.tetrahedra[4];
                pyraminx->tetrahedra[4] = temp.tetrahedra[7];
                pyraminx->tetrahedra[7] = temp.tetrahedra[1];
            }
        }
        if(axis == 3) {
            if(direction == 1) {
                pyraminx->tetrahedra[1] = temp.tetrahedra[3];
                pyraminx->tetrahedra[3] = temp.tetrahedra[6];
                pyraminx->tetrahedra[6] = temp.tetrahedra[1];
            }
            if(direction == -1) {
                pyraminx->tetrahedra[1] = temp.tetrahedra[6];
                pyraminx->tetrahedra[3] = temp.tetrahedra[1];
                pyraminx->tetrahedra[6] = temp.tetrahedra[3];
            }
        }
    }
    return;
}




Pyraminx rotatePyraminx(Pyraminx pyraminx, double turn_spd, int axis, int layer, int direction) {

    Tetrahedron **turning_tetrahedra = malloc(6 * sizeof(Tetrahedron *));
    Octahedron **turning_octahedra = malloc(3 * sizeof(Octahedron *));
    int total_tetrahedra;
    int total_octahedra;
    if(layer == 0) {
        total_tetrahedra = 6;
        total_octahedra = 3;
        if(axis == 0) {
            turning_tetrahedra[0] = &pyraminx.tetrahedra[0];
            turning_tetrahedra[1] = &pyraminx.tetrahedra[1];
            turning_tetrahedra[2] = &pyraminx.tetrahedra[2];
            turning_tetrahedra[3] = &pyraminx.tetrahedra[3];
            turning_tetrahedra[4] = &pyraminx.tetrahedra[4];
            turning_tetrahedra[5] = &pyraminx.tetrahedra[5];
            turning_octahedra[0] = &pyraminx.octahedra[0];
            turning_octahedra[1] = &pyraminx.octahedra[1];
            turning_octahedra[2] = &pyraminx.octahedra[2];
        }
        else if(axis == 1) {
            turning_tetrahedra[0] = &pyraminx.tetrahedra[0];
            turning_tetrahedra[1] = &pyraminx.tetrahedra[1];
            turning_tetrahedra[2] = &pyraminx.tetrahedra[2];
            turning_tetrahedra[3] = &pyraminx.tetrahedra[6];
            turning_tetrahedra[4] = &pyraminx.tetrahedra[7];
            turning_tetrahedra[5] = &pyraminx.tetrahedra[9];
            turning_octahedra[0] = &pyraminx.octahedra[0];
            turning_octahedra[1] = &pyraminx.octahedra[1];
            turning_octahedra[2] = &pyraminx.octahedra[3];
        }
        else if(axis == 2) {
            turning_tetrahedra[0] = &pyraminx.tetrahedra[0];
            turning_tetrahedra[1] = &pyraminx.tetrahedra[3];
            turning_tetrahedra[2] = &pyraminx.tetrahedra[5];
            turning_tetrahedra[3] = &pyraminx.tetrahedra[6];
            turning_tetrahedra[4] = &pyraminx.tetrahedra[8];
            turning_tetrahedra[5] = &pyraminx.tetrahedra[9];
            turning_octahedra[0] = &pyraminx.octahedra[0];
            turning_octahedra[1] = &pyraminx.octahedra[2];
            turning_octahedra[2] = &pyraminx.octahedra[3];
            
        }
        else if(axis == 3) {
            turning_tetrahedra[0] = &pyraminx.tetrahedra[2];
            turning_tetrahedra[1] = &pyraminx.tetrahedra[4];
            turning_tetrahedra[2] = &pyraminx.tetrahedra[5];
            turning_tetrahedra[3] = &pyraminx.tetrahedra[7];
            turning_tetrahedra[4] = &pyraminx.tetrahedra[8];
            turning_tetrahedra[5] = &pyraminx.tetrahedra[9];
            turning_octahedra[0] = &pyraminx.octahedra[1];
            turning_octahedra[1] = &pyraminx.octahedra[2];
            turning_octahedra[2] = &pyraminx.octahedra[3];
        }
    }
    else if(layer == 1) {
        total_tetrahedra = 3;
        total_octahedra = 1;
        if(axis == 0) {
            turning_tetrahedra[0] = &pyraminx.tetrahedra[6];
            turning_tetrahedra[1] = &pyraminx.tetrahedra[7];
            turning_tetrahedra[2] = &pyraminx.tetrahedra[8];
            turning_octahedra[0] = &pyraminx.octahedra[3];
        }
        else if(axis == 1) {
            turning_tetrahedra[0] = &pyraminx.tetrahedra[3];
            turning_tetrahedra[1] = &pyraminx.tetrahedra[4];
            turning_tetrahedra[2] = &pyraminx.tetrahedra[8];
            turning_octahedra[0] = &pyraminx.octahedra[2];
        }
        else if(axis == 2) {
            turning_tetrahedra[0] = &pyraminx.tetrahedra[1];
            turning_tetrahedra[1] = &pyraminx.tetrahedra[4];
            turning_tetrahedra[2] = &pyraminx.tetrahedra[7];
            turning_octahedra[0] = &pyraminx.octahedra[1];
            
        }
        else if(axis == 3) {
            turning_tetrahedra[0] = &pyraminx.tetrahedra[1];
            turning_tetrahedra[1] = &pyraminx.tetrahedra[3];
            turning_tetrahedra[2] = &pyraminx.tetrahedra[6];
            turning_octahedra[0] = &pyraminx.octahedra[0];
        }

    }
    else if(layer == 2) {
        total_tetrahedra = 1;
        total_octahedra = 0;
        if(axis == 0) {
            turning_tetrahedra[0] = &pyraminx.tetrahedra[9];
        }
        else if(axis == 1) {
            turning_tetrahedra[0] = &pyraminx.tetrahedra[5];
        }
        else if(axis == 2) {
            turning_tetrahedra[0] = &pyraminx.tetrahedra[2];
        }
        else if(axis == 3) {
            turning_tetrahedra[0] = &pyraminx.tetrahedra[0];
        }
    }
    turn_spd *= direction;
    pyraminx.layer_angle[axis][layer] += turn_spd;
    Matrix rotation = rotation_axis(pyraminx.axis_vectors[axis], turn_spd);

    for(int i = 0; i < total_tetrahedra; i++) {
        for(int v = 0; v < 4; v++) 
        {
            turning_tetrahedra[i]->vertices[v] = vertex_rotate(turning_tetrahedra[i]->vertices[v], rotation);
            *(turning_tetrahedra[i]) = updateTetrahedronEdgesAndFaces(*(turning_tetrahedra[i]));
        }
    }

    for(int i = 0; i < total_octahedra; i++) {
        for(int v = 0; v < 6; v++) 
        {
            turning_octahedra[i]->vertices[v] = vertex_rotate(turning_octahedra[i]->vertices[v], rotation);
            *(turning_octahedra[i]) = updateOctahedronEdgesAndFaces(*(turning_octahedra[i]));
        }
    }



    free(turning_tetrahedra);
    free(turning_octahedra);

    return pyraminx;
}



































































int main(int argc, char* argv[]) {


    Pixel **pixels = malloc((WINDOW_WIDTH) * sizeof(Pixel*));
    for(int i = 0; i < WINDOW_WIDTH; i++) {
        pixels[i] = malloc((WINDOW_HEIGHT) * sizeof(Pixel));
    }




    int pixel_size = 1;
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        printf("Could not initialize SDL: %s\n", SDL_GetError());
        return 1;
    }

    SDL_Window* window = SDL_CreateWindow("3D Renderer", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, WINDOW_WIDTH*pixel_size, WINDOW_HEIGHT*pixel_size, 0);
    if (!window) {
        printf("Could not create window: %s\n", SDL_GetError());
        SDL_Quit();
        return 1;
    }

    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_SOFTWARE);
    if (!renderer) {
        printf("Could not create renderer: %s\n", SDL_GetError());
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }

    SDL_Surface *surface = SDL_CreateRGBSurface(0, WINDOW_WIDTH, WINDOW_HEIGHT, 32, 0x00FF0000, 0x0000FF00, 0x000000FF, 0xFF000000);

    SDL_Event event;
    bool running = true;

    double speed_scale = 1.0;

    // Rubiks Cube
    RubiksCube rubiksCube;
    Vertex rubiksCubePosition = {-3.0, 2.5, 0.0};
    double rubiksCubeSize = 3.0 / (double)RUBIKS_CUBE_DIM;
    rubiksCube = createRubiksCube(rubiksCube, rubiksCubeSize);
    double turn_spd = degrees_to_radians(5);
    double turn_snap_bias = turn_spd - 0.1;

    // Player
    Vertex playerPosition = {-6.0, 6.0, 6.0};
    double playerYSpd = 0.5;
    double playerHSpd = 0.5;

    double yaw = degrees_to_radians(-45);
    double pitch = degrees_to_radians(-30);
    double camera_spd = degrees_to_radians(5);
    double camera_snap_bias = camera_spd - 0.1; // bias for snap to perpendicular angle

    // light source
    Vertex lightSourcePosition = {-5, 5, -5};






    // floor grid
    FloorGrid floorGrid;
    floorGrid = createFloorGrid(floorGrid, 1.0, FLOOR_GRID_COLOR);
    Square floorGridSquare;
    floorGridSquare = createSquare((Vertex){0, -.05, 0}, FLOOR_GRID_DIM, FLOOR_GRID_COLOR);//-2.77
    floorGridObjPtr = malloc(sizeof(uint32_t));
    floorGridObjPtr = (uint32_t *)(&floorGridSquare);



    // tetrahedron
    double edgeLength = 1.0;
    Vertex tetrahedronPosition = {0.0, 0.0, 0.0};
    Tetrahedron tetrahedron = createTetrahedron(tetrahedronPosition, edgeLength);
    
    
    
    
    // pyraminx
    Pyraminx pyraminx;
    pyraminx = createPyraminx(pyraminx);
    double pyraminx_turn_spd = degrees_to_radians(5);
    double pyraminx_turn_snap_bias = pyraminx_turn_spd - 0.1;
    int temp_axis = 0;
    Vertex pyraminxPosition = {3.0, 1.75, 0.0};









    // collision array
    collisionTris.tri_address = malloc(MAX_TRIS * sizeof(Vertex *));
    collisionTris.is_quad_tri = malloc(MAX_TRIS * sizeof(bool));
    collisionTris.object_address = malloc(MAX_TRIS * sizeof(uint32_t *));

    collisionTris.v = malloc(MAX_TRIS * sizeof(Vertex *));
    for(int i = 0; i < MAX_TRIS; i++) {
        collisionTris.v[i] = malloc(3 * sizeof(Vertex));
    }

    




    // key stuff
    bool key_down = false;
    int key_count = 10;

    // rubiks cube turning stuff
    bool turn_cube = false;
    int axis = -1;
    int layer = -1;
    int direction = 0;

    // pyraminx turning stuff
    bool turn_pyraminx = false;


    int mode = 0;
    bool draw_floor_grid = true;

    // loop
    while(running == true) 
    {

    int TARGET_FPS;
    if(drawShadows == false) {
        TARGET_FPS = 20;
    }
    else if(drawShadows = true) {
        TARGET_FPS = 10;
    }
    int TARGET_FRAME_DURATION = 1000 / TARGET_FPS;


        Uint64 frameStart = SDL_GetTicks64();
        while (SDL_PollEvent(&event) != 0) {
            if (event.type == SDL_QUIT) {
                running = false;
            }
        }
        const Uint8* key_state = SDL_GetKeyboardState(&key_count);

        // quit if the user presses escape key
        if (key_state[SDL_SCANCODE_ESCAPE]) {
            running = false;
            continue;
        }




        // change puzzle mode
        if (key_state[SDL_SCANCODE_Q]) {
            mode = 0; // rubiks cube mode
            printf("MODE: CUBE\n");
        }
        else if (key_state[SDL_SCANCODE_P]) {
            mode = 1; // pyraminx mode
            printf("MODE: PYRAMINX\n");
        }

        // change puzzle turning axis
        if(key_state[SDL_SCANCODE_1]) {
            temp_axis = 0;
            printf("AXIS: 1\n");
        }
        else if(key_state[SDL_SCANCODE_2]) {
            temp_axis = 1;
            printf("AXIS: 2\n");
        }
        else if(key_state[SDL_SCANCODE_3]) {
            temp_axis = 2;
            printf("AXIS: 3\n");
        }
        else if(key_state[SDL_SCANCODE_4] && mode == 1) {
            temp_axis = 3;
            printf("AXIS: 4\n");
        }

        // rubiks cube turning inputs
        if(turn_cube == false && mode == 0) {
            if (key_state[SDL_SCANCODE_B]) {
                axis = temp_axis;
                layer = 0;
                direction = 1;
                turn_cube = true;
            }
            else if (key_state[SDL_SCANCODE_N]) {
                axis = temp_axis;
                layer = 0;
                direction = -1;
                turn_cube = true;
            }
            else if (key_state[SDL_SCANCODE_G]) {
                axis = temp_axis;
                layer = 1;
                direction = 1;
                turn_cube = true;
            }
            else if (key_state[SDL_SCANCODE_H]) {
                axis = temp_axis;
                layer = 1;
                direction = -1;
                turn_cube = true;
            }
            else if (key_state[SDL_SCANCODE_T]) {
                axis = temp_axis;
                layer = 2;
                direction = 1;
                turn_cube = true;
            }
            else if (key_state[SDL_SCANCODE_Y]) {
                axis = temp_axis;
                layer = 2;
                direction = -1;
                turn_cube = true;
            }
            else if (key_state[SDL_SCANCODE_5]) {
                axis = temp_axis;
                layer = 3;
                direction = 1;
                turn_cube = true;
            }
            else if (key_state[SDL_SCANCODE_6]) {
                axis = temp_axis;
                layer = 3;
                direction = -1;
                turn_cube = true;
            }
        }



        // pyraminx turning inputs
        if(turn_pyraminx == false && mode == 1) {
            
            if (key_state[SDL_SCANCODE_B]) {
                axis = temp_axis;
                layer = 0;
                direction = 1;
                turn_pyraminx = true;
            }
            else if (key_state[SDL_SCANCODE_N]) {
                axis = temp_axis;
                layer = 0;
                direction = -1;
                turn_pyraminx = true;
            }
            else if (key_state[SDL_SCANCODE_G]) {
                axis = temp_axis;
                layer = 1;
                direction = 1;
                turn_pyraminx = true;
            }
            else if (key_state[SDL_SCANCODE_H]) {
                axis = temp_axis;
                layer = 1;
                direction = -1;
                turn_pyraminx = true;
            }
            else if (key_state[SDL_SCANCODE_T]) {
                axis = temp_axis;
                layer = 2;
                direction = 1;
                turn_pyraminx = true;
            }
            else if (key_state[SDL_SCANCODE_Y]) {
                axis = temp_axis;
                layer = 2;
                direction = -1;
                turn_pyraminx = true;
            }
        }









        






        // rotate camera
        if (key_state[SDL_SCANCODE_LEFT]) {
            yaw += camera_spd;
            yaw = fix_angle(yaw, camera_snap_bias);
        }
        if (key_state[SDL_SCANCODE_RIGHT]) {
            yaw -= camera_spd;
            yaw = fix_angle(yaw, camera_snap_bias);
        }
        if (key_state[SDL_SCANCODE_UP]) {
            pitch += camera_spd;
            pitch = fix_angle(pitch, camera_snap_bias);
        }
        if (key_state[SDL_SCANCODE_DOWN]) {
            pitch -= camera_spd;
            pitch = fix_angle(pitch, camera_snap_bias);
        }

















        double player_h_spd = 0.5;
        double player_y_spd = 0.5;
        if ((key_state[SDL_SCANCODE_W] && key_state[SDL_SCANCODE_A]) ^ (key_state[SDL_SCANCODE_A] && key_state[SDL_SCANCODE_S]) ^ (key_state[SDL_SCANCODE_S] && key_state[SDL_SCANCODE_D]) ^ (key_state[SDL_SCANCODE_D] && key_state[SDL_SCANCODE_W])) {
           player_h_spd /= sqrt(2);
        }
        double sinYaw = sin(yaw);
        double cosYaw = cos(yaw);
        double sinPitch = sin(pitch);
        double cosPitch = cos(pitch);

        // move player
        if (key_state[SDL_SCANCODE_W]) {
            playerPosition.x -= player_h_spd * sinYaw;
            playerPosition.z -= player_h_spd * cosYaw;
        }
        if (key_state[SDL_SCANCODE_A]) {
            playerPosition.x -= player_h_spd * cosYaw;
            playerPosition.z += player_h_spd * sinYaw;
        }
        if (key_state[SDL_SCANCODE_S]) {
            playerPosition.x += player_h_spd * sinYaw;
            playerPosition.z += player_h_spd * cosYaw;
        }
        if (key_state[SDL_SCANCODE_D]) {
            playerPosition.x += player_h_spd * cosYaw;
            playerPosition.z -= player_h_spd * sinYaw;
        }
        if (key_state[SDL_SCANCODE_SPACE]) {
            playerPosition.y += player_y_spd;
        }
        if (key_state[SDL_SCANCODE_TAB]) {
            playerPosition.y -= player_y_spd;
        }








        double light_source_h_spd = 0.3;
        double light_source_y_spd = 0.5;
        if ((key_state[SDL_SCANCODE_8] && key_state[SDL_SCANCODE_4]) ^ (key_state[SDL_SCANCODE_4] && key_state[SDL_SCANCODE_5]) ^ (key_state[SDL_SCANCODE_5] && key_state[SDL_SCANCODE_6]) ^ (key_state[SDL_SCANCODE_6] && key_state[SDL_SCANCODE_8])) {
            light_source_h_spd /= sqrt(2);
        }



        // move light source
        if (key_state[SDL_SCANCODE_EQUALS]) {
            lightSourcePosition.x -= light_source_h_spd * sinYaw;
            lightSourcePosition.z -= light_source_h_spd * cosYaw;
        }
        if (key_state[SDL_SCANCODE_LEFTBRACKET]) {
            lightSourcePosition.x -= light_source_h_spd * cosYaw;
            lightSourcePosition.z += light_source_h_spd * sinYaw;
        }
        if (key_state[SDL_SCANCODE_RIGHTBRACKET]) {
            lightSourcePosition.x += light_source_h_spd * sinYaw;
            lightSourcePosition.z += light_source_h_spd * cosYaw;
        }
        if (key_state[SDL_SCANCODE_BACKSLASH]) {
            lightSourcePosition.x += light_source_h_spd * cosYaw;
            lightSourcePosition.z -= light_source_h_spd * sinYaw;
        }
        if (key_state[SDL_SCANCODE_BACKSPACE]) {
            lightSourcePosition.y += light_source_y_spd;
        }
        if (key_state[SDL_SCANCODE_RETURN]) {
            lightSourcePosition.y -= light_source_y_spd;
        }



        // turn shadows on/off
        if (key_state[SDL_SCANCODE_PERIOD]) {
            drawShadows = false;
        }
        else if (key_state[SDL_SCANCODE_SLASH]) {
            drawShadows = true;
        }







        Matrix rotY = rotation_y(yaw);
        Matrix rotX = rotation_x(pitch);
        Matrix invRotY = rotation_y(-yaw);
        Matrix invRotX = rotation_x(-pitch);

        double bias = 0.01;
        if (mode == 0 && turn_cube == true) {
            rubiksCube = rotateRubiksCube(rubiksCube, turn_spd, axis, layer, direction);
            if (is_fixed_90(rubiksCube.layer_angle[axis][layer], bias)) {
                updateCubiesIndexing(&rubiksCube, axis, layer, direction);
                int axis = -1;
                int layer = -1;
                int direction = 0;
                turn_cube = false;
            }
        }
        else if (mode == 1 && turn_pyraminx == true) {
            bias = 0.05;
            pyraminx = rotatePyraminx(pyraminx, pyraminx_turn_spd, axis, layer, direction);
            if (is_fixed_120(pyraminx.layer_angle[axis][layer], bias)) {
                updatePyraminxIndexing(&pyraminx, axis, layer, direction);
                int axis = -1;
                int layer = -1;
                int direction = 0;
                turn_pyraminx = false;
            }
        }




        if(drawShadows == true || drawShadows == false) {
            collision_n = 0;
            for(int i = 0; i < RUBIKS_CUBE_DIM; i++) {
                for(int j = 0; j < RUBIKS_CUBE_DIM; j++) {
                    for(int k = 0; k < RUBIKS_CUBE_DIM; k++) {
                        for(int f = 0; f < 12; f++) {
                            bool is_correct_layer = false;
                            if(turn_cube == true) {
                                switch(axis) {
                                    case 0:
                                        is_correct_layer = (i == layer);
                                        break;
                                    case 1:
                                        is_correct_layer = (j == layer);
                                        break;
                                    case 2:
                                        is_correct_layer = (k == layer);
                                        break;
                                    default:
                                        is_correct_layer = false;
                                        break;
                                }
                            }

                            bool is_face_in_turning_layer = isFaceInTurningLayer(rubiksCube.cubies[i][j][k].faces[f].tri, axis, layer, rubiksCubePosition, rubiksCubeSize);

                            if(isOuterFace(rubiksCube.cubies[i][j][k].faces[f].tri, rubiksCubeSize, turn_cube, is_correct_layer, is_face_in_turning_layer, axis, layer, rubiksCube, (int *)NULL) == true) {
                                collisionTris.tri_address[collision_n] = (Vertex *)rubiksCube.cubies[i][j][k].faces[f].tri;

                                if(f % 2 == 0) {
                                    // even
                                    collisionTris.is_quad_tri[collision_n] = true;
                                }
                                else {
                                    // odd
                                    collisionTris.is_quad_tri[collision_n] = true;
                                }

                                collisionTris.object_address[collision_n] = (uint32_t *)(&rubiksCube);

                                collisionTris.v[collision_n][0] = vertex_add(rubiksCube.cubies[i][j][k].faces[f].tri[0], rubiksCubePosition);
                                collisionTris.v[collision_n][1] = vertex_add(rubiksCube.cubies[i][j][k].faces[f].tri[1], rubiksCubePosition);
                                collisionTris.v[collision_n][2] = vertex_add(rubiksCube.cubies[i][j][k].faces[f].tri[2], rubiksCubePosition);
                                collision_n++;
                            }   
                        }
                    }
                }
            }

            for(int i = 0; i < TOTAL_TETRAHEDRA; i++) {
                for(int f = 0; f < 4; f++) {
                    bool is_correct_layer = false;
                    if(turn_pyraminx == true) 
                    {
                        if(axis == 0) 
                        {
                            if(layer == 0) {
                                if(i == 0 || i == 1 || i == 2 || i == 3 || i == 4 || i == 5) {
                                    is_correct_layer = true;
                                }
                            }
                            else if(layer == 1) {
                                if(i == 6 || i == 7 || i == 8) {
                                    is_correct_layer = true;
                                }
                            }
                            else if(layer == 2) {
                                if(i == 9) {
                                    is_correct_layer = true;
                                }
                            }
                        }

                        else if(axis == 1) 
                        {
                            if(layer == 0) {
                                if(i == 0 || i == 1 || i == 2 || i == 6 || i == 7 || i == 9) {
                                    is_correct_layer = true;
                                }
                            }
                            else if(layer == 1) {
                                if(i == 3 || i == 4 || i == 8) {
                                    is_correct_layer = true;
                                }
                            }
                            else if(layer == 2) {
                                if(i == 5) {
                                    is_correct_layer = true;
                                }
                            }
                        }
                        else if(axis == 2) 
                        {
                            if(layer == 0) {
                                if(i == 0 || i == 3 || i == 5 || i == 6 || i == 8 || i == 9) {
                                    is_correct_layer = true;
                                }
                            }
                            else if(layer == 1) {
                                if(i == 1 || i == 4 || i == 7) {
                                    is_correct_layer = true;
                                }
                            }
                            else if(layer == 2) {
                                if(i == 2) {
                                    is_correct_layer = true;
                                }
                            }
                        }

                        else if(axis == 3) 
                        {
                            if(layer == 0) {
                                if(i == 2 || i == 4 || i == 5 || i == 7 || i == 8 || i == 9) {
                                    is_correct_layer = true;
                                }
                            }
                            else if(layer == 1) {
                                if(i == 1 || i == 3 || i == 6) {
                                    is_correct_layer = true;
                                }
                            }
                            else if(layer == 2) {
                                if(i == 0) {
                                    is_correct_layer = true;
                                }
                            }
                        }
                    }
                    if(isOuterPyraminxFace(pyraminx.tetrahedra[i].faces[f].tri, turn_pyraminx, is_correct_layer, axis, layer, pyraminx, (int *)NULL)) {
                        collisionTris.tri_address[collision_n] = (Vertex *)pyraminx.tetrahedra[i].faces[f].tri;
                        collisionTris.is_quad_tri[collision_n] = false;
                        collisionTris.object_address[collision_n] = (uint32_t *)(&pyraminx);
                        collisionTris.v[collision_n][0] = vertex_add(pyraminx.tetrahedra[i].faces[f].tri[0], pyraminxPosition);
                        collisionTris.v[collision_n][1] = vertex_add(pyraminx.tetrahedra[i].faces[f].tri[1], pyraminxPosition);
                        collisionTris.v[collision_n][2] = vertex_add(pyraminx.tetrahedra[i].faces[f].tri[2], pyraminxPosition);
                        collision_n++;
                        //collision[collision_n] = pyraminx.tetrahedra[i].faces[f].tri;
                        /*Vertex temp = collision[collision_n][1];
                        collision[collision_n][1] = collision[collision_n][2];
                        collision[collision_n][2] = temp;*/
                        //collision_n++;
                    }
                }
            }

            for(int i = 0; i < TOTAL_OCTAHEDRA; i++) {
                for(int f = 0; f < 8; f++) {
                    bool is_correct_layer = false;
                    if(turn_pyraminx == true) 
                    {
                        if(axis == 0) 
                        {
                            if(layer == 0) {
                                if(i == 0 || i == 1 || i == 2) {
                                    is_correct_layer = true;
                                }
                            }
                            else if(layer == 1) {
                                if(i == 3) {
                                    is_correct_layer = true;
                                }
                            }
                        }

                        else if(axis == 1) 
                        {
                            if(layer == 0) {
                                if(i == 0 || i == 1 || i == 3) {
                                    is_correct_layer = true;
                                }
                            }
                            else if(layer == 1) {
                                if(i == 2) {
                                    is_correct_layer = true;
                                }
                            }
                        }

                        else if(axis == 2) 
                        {
                            if(layer == 0) {
                                if(i == 0 || i == 2 || i == 3) {
                                    is_correct_layer = true;
                                }
                            }
                            else if(layer == 1) {
                                if(i == 1) {
                                    is_correct_layer = true;
                                }
                            }
                        }

                        else if(axis == 3) 
                        {
                            if(layer == 0) {
                                if(i == 1 || i == 2 || i == 3) {
                                    is_correct_layer = true;
                                }
                            }
                            else if(layer == 1) {
                                if(i == 0) {
                                    is_correct_layer = true;
                                }
                            }
                        }
                    }
                    if(isOuterPyraminxFace(pyraminx.octahedra[i].faces[f].tri, turn_pyraminx, is_correct_layer, axis, layer, pyraminx, (int *)NULL)) {
                        collisionTris.tri_address[collision_n] = (Vertex *)pyraminx.octahedra[i].faces[f].tri;
                        collisionTris.is_quad_tri[collision_n] = false;
                        collisionTris.object_address[collision_n] = (uint32_t *)(&pyraminx);
                        collisionTris.v[collision_n][0] = vertex_add(pyraminx.octahedra[i].faces[f].tri[0], pyraminxPosition);
                        collisionTris.v[collision_n][1] = vertex_add(pyraminx.octahedra[i].faces[f].tri[2], pyraminxPosition);
                        collisionTris.v[collision_n][2] = vertex_add(pyraminx.octahedra[i].faces[f].tri[1], pyraminxPosition);
                        collision_n++;
                    }
                }
            }

            if(draw_floor_grid == true) {
                for(int f = 0; f < 2; f++) {
                    collisionTris.tri_address[collision_n] = (Vertex *)floorGridSquare.faces[f].tri;
                    if(f % 2 == 0) {
                        collisionTris.is_quad_tri[collision_n] = true;
                    }
                    else {
                        collisionTris.is_quad_tri[collision_n] = true;
                    }
                    collisionTris.object_address[collision_n] = (uint32_t *)(&floorGridSquare);
                    collisionTris.v[collision_n][0] = floorGridSquare.faces[f].tri[0];
                    collisionTris.v[collision_n][1] = floorGridSquare.faces[f].tri[1];
                    collisionTris.v[collision_n][2] = floorGridSquare.faces[f].tri[2];
                    collision_n++;
                }
            }
        }

        




        // clear `pixels`
        for(int i = 0; i < WINDOW_WIDTH; i++) {
            for(int j = 0; j < WINDOW_HEIGHT; j++) {
                pixels[i][j].color = BACKGROUND_COLOR_DARK;
                pixels[i][j].depth = 100.0;
            }
        }

        // rubiks cube tris
        drawRubiksCubeTris(pixels, &rubiksCube, rubiksCubePosition, rubiksCubeSize, playerPosition, lightSourcePosition, rotX, rotY, invRotX, invRotY, turn_cube, axis, layer);

        // pyraminx tris
        drawPyraminxTris(pixels, &pyraminx, pyraminxPosition, playerPosition, lightSourcePosition, rotX, rotY, invRotX, invRotY, turn_pyraminx, axis, layer);

        // remove/add floor grid
        if(key_state[SDL_SCANCODE_0]) { 
            draw_floor_grid = 1;
        }
        if(key_state[SDL_SCANCODE_9]) { 
            draw_floor_grid = 0;
        }
        // floor grid tris
        if(draw_floor_grid) {
           drawFloorGridTris(pixels, &floorGridSquare, playerPosition, lightSourcePosition, rotX, rotY, invRotX, invRotY);
        }


        // rubiks cube lines
        drawRubiksCubeLines(pixels, &rubiksCube, rubiksCubePosition, rubiksCubeSize, playerPosition, lightSourcePosition, rotX, rotY, invRotX, invRotY, turn_cube, axis, layer);

        // pyraminx lines
        drawPyraminxLines(pixels, &pyraminx, pyraminxPosition, playerPosition, lightSourcePosition, rotX, rotY, invRotX, invRotY, turn_pyraminx, axis, layer);

        // floor grid lines
        if(draw_floor_grid) {
            drawFloorGridLines(pixels, &floorGrid, &floorGridSquare, playerPosition, lightSourcePosition, rotX, rotY, invRotX, invRotY);
        }

        // light source cube
        double cube_size = 0.5;
        drawLightSourceCubeTris(pixels, cube_size, playerPosition, lightSourcePosition, rotX, rotY, invRotX, invRotY);
        drawLightSourceCubeLines(pixels, cube_size, playerPosition, lightSourcePosition, rotX, rotY, invRotX, invRotY);
        

        // copy pixel data to the SDL_Surface
        int *surfacePixels = (int *)surface->pixels;
        for (int y = 0; y < WINDOW_HEIGHT; ++y) {
            for (int x = 0; x < WINDOW_WIDTH; ++x) {
                surfacePixels[y * WINDOW_WIDTH + x] = pixels[x][y].color;
            }
        }

        SDL_Texture *texture = SDL_CreateTextureFromSurface(renderer, surface);
        SDL_RenderClear(renderer);
        SDL_RenderCopy(renderer, texture, NULL, NULL);
        SDL_RenderPresent(renderer);

        SDL_DestroyTexture(texture);

        Uint64 frameEnd = SDL_GetTicks64();
        Uint64 frameDuration = frameEnd - frameStart;
        if (frameDuration < TARGET_FRAME_DURATION) {
            SDL_Delay(TARGET_FRAME_DURATION - frameDuration);
        }

        Uint64 newFrameEnd = SDL_GetTicks64();
        Uint64 newFrameDuration = newFrameEnd - frameStart;
        double fps = 1000.0 / newFrameDuration;
        //printf("FPS: %d\n", (int)fps);
    }

    for(int i = 0; i < WINDOW_WIDTH; i++) {
        free(pixels[i]);
    }
    free(pixels);


    for(int i = 0; i < MAX_TRIS; i++) {
        free(collisionTris.v[i]);
    }

    free(collisionTris.object_address);
    free(collisionTris.is_quad_tri);
    free(collisionTris.v);
    free(collisionTris.tri_address);

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}
