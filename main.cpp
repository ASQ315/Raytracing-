#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <utility>

using namespace std;

class Vectores {
    public:
        double x, y, z;
        Vectores () {};
        Vectores (double x, double y, double z) : x(x), y(y), z(z) {}

        Vectores operator+ (const Vectores& v) const {
            return Vectores(x + v.x, y + v.y, z + v.z);
        }

        Vectores operator- (const Vectores& v) const {
            return Vectores(x - v.x, y - v.y, z - v.z);
        }

        Vectores operator* (double s) const {
            return Vectores(x * s, y * s, z * s);
        }

        double dot (const Vectores& v) const {
            return x * v.x + y * v.y + z * v.z;
        }
        
        Vectores normalize() const {
            double mag = sqrt(x * x + y * y + z * z);
            return Vectores(x / mag, y / mag, z / mag);
        }

};

class Objects {
    public:

        virtual ~Objects() = default;
};

class Sphere : public Objects{
    public:
        Vectores Center;
        double Radius;
        Vectores Color;
        int Specular;
        double Reflective;

        Sphere(Vectores center, double radius, Vectores color, int specular, double reflective) : 
            Center(center), Radius(radius), Color(color), Specular(specular), Reflective(reflective) {}
               
};

class Plane : public Objects{
    public:
        Vectores Point;
        Vectores Normal;
        Vectores Color;
        int Specular;
        double Reflective;

        Plane(Vectores point, Vectores normal, Vectores color, int specular, double reflective) :
            Point(point), Normal(normal), Color(color), Specular(specular), Reflective(reflective) {}
            
        };

class Light : public Objects {
    public:
        string Type;
        double Intensity;
        Vectores Position;
        Vectores Direction;

        Light(string type, double intensity, Vectores position, Vectores direction) :
            Type(type), Intensity(intensity), Position(position), Direction(direction) {}
            
        };

vector<Sphere> spheres;
vector<Light> lights;
vector<Plane> planes;

int clamp(double value, double min_val, double max_val) {
    return static_cast <int> (max(min_val, min(max_val, value)));
}

double length(const Vectores& v) {
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

bool IntersectRaySphere(const Vectores& O, const Vectores& D, const Sphere& sphere, double& t1, double& t2) {
    Vectores CO = O - sphere.Center;

    double a = D.dot(D);
    double b = 2 * CO.dot(D);
    double c = CO.dot(CO) - sphere.Radius * sphere.Radius;

    double discriminant = b * b - 4 * a * c;

    if (discriminant < 0) {
        return false;
    }

    t1 = (-b + sqrt(discriminant)) / (2 * a);
    t2 = (-b - sqrt(discriminant)) / (2 * a);

    return true;
}

double IntersectRayPlane(const Vectores& O,const Vectores&  D,const Plane& plane) {

    Vectores normal = plane.Normal;
    double denom = normal.dot(D);

    if (abs(denom) > 1e-6 ){
        double t = (plane.Point - O).dot(normal) / denom;
        
        if (t > 0) {
            return t; 
        }
        else {
            return INFINITY;
        }
    }
    return INFINITY;
}

pair<const Objects*, double> ClosestIntersection(const Vectores& O, const Vectores& D, double t_min, double t_max) {
    double closest_t = INFINITY;
    const Objects* closest_object = nullptr;

    for (const Sphere& sphere : spheres) {
        double t1, t2;

        if (IntersectRaySphere(O, D, sphere, t1, t2)) {
            if(t1 >= t_min && t1 <= t_max && t1 < closest_t) {
                closest_t = t1;
                closest_object = &sphere;
            }
            if(t2 >= t_min && t2 <= t_max && t2 < closest_t) {
                closest_t = t2;
                closest_object = &sphere;
            }
        }
    }

    for (const Plane& plane : planes) {
        double t = IntersectRayPlane(O, D, plane);

        if(t_min <= t && t <= t_max && t < closest_t) {
            closest_t = t;
            closest_object = &plane;
        }
    }

    return {closest_object, closest_t};
}

Vectores ReflectRay(const Vectores& R, const Vectores& N) {
    return  N * 2.0 * N.dot(R) - R;
}

double ComputeLighting(const Vectores& P, const Vectores& N, const Vectores& V, const int s) {
    double i = 0.0;
    Vectores L;
    double t_max;

    for (const Light& light : lights) {
        if (light.Type == "ambient") {
          i += light.Intensity; 
        }
        else {
            if(light.Type == "point") {
               L =  (light.Position - P).normalize();
               t_max = length(light.Position - P);

            }
            else {
               L = light.Direction.normalize(); 
               t_max = INFINITY;

            }
            //Shadow Check
            pair<const Objects*, double> shadow_result = ClosestIntersection(P, L, 0.001, t_max);

            if (shadow_result.first != nullptr) {
                continue;
            }
            //Diffuse
            double n_dot_l = N.dot(L);
            if(n_dot_l > 0) {
                i += light.Intensity * n_dot_l;
            }
            //Specular
            if (s != -1) {
                Vectores R = ReflectRay(L, N);
                double r_dot_v = R.dot(V);

                if( r_dot_v > 0) {
                    i += light.Intensity * pow(r_dot_v / (length(R) * length(V)), s);
                }
            }
        }

    }
    return i;
}

Vectores TraceRay(const Vectores O, const Vectores D, double t_min, double t_max, int recursion_depth){
    
    pair <const Objects*, double> intersection_result = ClosestIntersection(O, D, t_min, t_max);
    const Objects* closest_object = intersection_result.first;
    double closest_t = intersection_result.second;

    if(closest_object == nullptr) {
        return Vectores(255, 255, 255); //rgb(255, 255, 255)Background Color
    }

    Vectores P = O + D * closest_t;
    Vectores N;
    Vectores object_color;
    int specular = -1;
    double reflective;
    
     if(const Sphere* sphere = dynamic_cast<const Sphere*> (closest_object)) {

        N = (P - sphere -> Center).normalize();
        object_color = sphere->Color;
        specular = sphere->Specular;
        reflective = sphere->Reflective;
        
    
    }
    else if(const Plane* plane = dynamic_cast<const Plane*> (closest_object)) {

        N = plane->Normal;
        object_color = plane->Color;
        specular = plane->Specular;
        reflective = plane->Reflective;

    }

    Vectores V = D * -1.0;
    double intensity = ComputeLighting(P, N, V, specular);
    Vectores local_color = object_color * intensity;
    
    if (recursion_depth <= 0 || reflective <= 0 ){
        return local_color;
    }

    Vectores R = ReflectRay(V, N);
    Vectores reflected_color = TraceRay(P, R, 0.001, INFINITY, recursion_depth - 1);
    
    Vectores final_color = local_color * (1 - reflective) + reflected_color * reflective;

    final_color.x = clamp(final_color.x, 0.0, 255.0);
    final_color.y = clamp(final_color.y, 0.0, 255.0);
    final_color.z = clamp(final_color.z, 0.0, 255.0);

    return final_color;
}
Vectores CanvasToViewport(int x, int y, int Cw, int Ch, int Vw, int Vh, int d) {
    return Vectores(
        (x + 0.5) * static_cast<double>(Vw) / Cw,
        (y + 0.5) * static_cast<double>(Vh) / Ch,
        d
    );
}

double smoothStep(double t) {
    return t * t * (3 - 2 * t); 
}

int main(){
    int Cw = 500;
    int Ch = 500;
    int Vw = 1;
    int Vh = 1;
    int d = 1;

    int recursion_depth = 3;

    Vectores O = Vectores(0, 0, 0);

    spheres = {
        Sphere(Vectores(0, 2, 7), 0.5, Vectores(31, 248, 103), 500, 0.7), //rgb(31, 248, 103)
        Sphere(Vectores(2, 0, 7), 0.5, Vectores(249, 254, 84), 500, 0.7), //rgb(249, 254, 84)
        Sphere(Vectores(-2, 0, 7), 0.5, Vectores(61, 55, 240), 500, 0.7), //rgb(61, 55, 240)
        Sphere(Vectores(0, -2, 7), 0.5, Vectores(243, 136, 230), 500, 0.7), //rgb(243, 136, 230)
    };

    planes = {
        Plane(Vectores(0, -3, 0), Vectores(0, 1, 0), Vectores(115, 82, 233), 500, 1),  // rgb(115, 82, 233)
        Plane(Vectores(0, 3, 0), Vectores(0, -1, 0), Vectores(115, 82, 233), 500, 1), 
        Plane(Vectores(-4, 3, 0), Vectores(-1, 0, 0), Vectores(8, 4, 70), 500, 1), // rgb(8, 4, 70) 
        Plane(Vectores(4, 3, 0), Vectores(1, 0, 0), Vectores(8, 4, 70), 500, 1)  
    }; 

    lights = {
        Light("ambient", 0.2, Vectores(), Vectores()),
        Light("point", 0.6, Vectores(0, 0, 0), Vectores()),
        Light("directional", 0.2, Vectores(), Vectores(1, 4, 4))
    };

    int num_frames = 120;
    int merge_phase = 50; 
    int swapped_phase = 80; 

Vectores start_positions[] = {
    Vectores(0, 2, 7),  
    Vectores(2, 0, 7),  
    Vectores(-2, 0, 7),
    Vectores(0, -2, 7)   
};

Vectores center_position(0, 0, 7); 

Vectores new_positions[] = {
    start_positions[3], 
    start_positions[1], 
    start_positions[2], 
    start_positions[0]  
};

for (int frame = 0; frame < num_frames; frame++) {
    double t;
    if (frame < merge_phase) {

        t = smoothStep((double)frame / (merge_phase - 1));
        for (int i = 0; i < 4; i++) {
            spheres[i] = Sphere(
                start_positions[i] * (1 - t) + center_position * t,
                0.5, spheres[i].Color, 500, 0
            );
        }
    }
    else {
        double phase_t = (double)(frame - merge_phase) / (swapped_phase - merge_phase);
        t = smoothStep(phase_t);
       
        for (int i = 0; i < 4; i++) {
            spheres[i] = Sphere(
                (center_position) * (1 - t) + new_positions[i] * t,
                0.5, spheres[i].Color, 500, 0
            );
        }
    }
        
        string filename = "frames/frame_" + to_string(frame) + ".ppm";
        ofstream filePPM(filename);
        filePPM << "P3\n" << Cw << " " << Ch << "\n255\n";

        for (int y = Ch / 2 - 1; y >= -Ch / 2; --y) {
            for (int x = -Cw / 2; x < Cw / 2; ++x) {
                Vectores D = CanvasToViewport(x, y, Cw, Ch, Vw, Vh, d);
                Vectores color = TraceRay(O, D, 1, INFINITY, recursion_depth);

                int r = clamp(color.x, 0.0, 255.0);
                int g = clamp(color.y, 0.0, 255.0);
                int b = clamp(color.z, 0.0, 255.0);

                filePPM << r << " " << g << " " << b << " ";
            }
            filePPM << "\n";
        }
        cout << "Frame " << frame << " generado" << endl;
    }
    return 0;
}