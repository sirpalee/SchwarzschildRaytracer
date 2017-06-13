#include <random>
#include <array>
#include <algorithm>
#include <numeric>
#include <future>
#include <mutex>
#include <type_traits>
#include <SFML/Graphics.hpp>
#include <SFML/OpenGL.hpp>

#include "vec.h"
#include "color.h"

static constexpr int diskWidth = 1024 + 512;
static constexpr int diskHeight = 256 * 8;

static const int nths = std::thread::hardware_concurrency();

struct Params {
    double M;
    double G;
    double fovw;//degrees
    double zcam;
    double ycam;
    double T0disk;//Kelvin
};

template<typename T>
struct ThetaPhi {
    T theta, phi;
};

struct Star {
    ThetaPhi<float> angles;
    Vec<float> colxyz;
    float d, T;
};

const int blocksize = 16;

struct TraceParams {
    int bx, by;
    int w, h;
    bool* finish;

    void reset() {
        //const auto q = h/blocksize;
        by = 0;//w/(blocksize*3);
        bx = 0;
    }

    void step() {
        if (!*finish) {
            const auto p = w / blocksize;
            const auto q = h / blocksize;

            if (by > q) {
                by = 0;
                *finish = true;
            }
            if (bx > p) {
                bx = 0;
                by += 1;
            }
            else { bx += 1; }
        }
    }

    float ratio_complete() const {
        const auto p = w / blocksize;
        const auto q = h / blocksize;
        return (by * p + bx) / ((float) (p + 2) * q);
    }
};

void trace_n(TraceParams& tpar, std::vector<Star> const& stars, std::vector<std::vector<size_t>> const& starmap,
             std::vector<Vec<float>>& xyz, std::vector<sf::Color>& res, std::vector<float> const& disk,
             Params const& params) {
    std::random_device rd;
    std::mutex mxyz;

    const double rs = 2.; //Schwarzschild radius
    const double robs = length(Vec<double>{0.0, params.ycam, params.zcam}); //Distance of observer
    const double rho_obs = sqrt(1. - rs / robs);
    const double da = params.fovw / (double) tpar.w;

    const auto rdin = 3.001 * rs; //Inner radius of disk in Schwarzschild units
    const auto rdout = 8.000 * rs; //Outer radius of disk in Schwarzschild units
    const auto hdisk = 0.040 * rs; //Disk thickness in Schwarzschild units

    const auto up = Vec<double>{0.0, 1.0, 0.0};

    const Vec<double> bh0 = {0.0, 0.0, 0.0};//{+1.5, 0.0, +3.15};
    // const Vec<double> bh1 = {-1.5, 0.0, -3.15};

    //Disk velocity vector at location
    auto vel_at = [&](Vec<double> const& p) {
        auto b1 = p - bh0;
        // auto b2 = p - bh1;
        auto R1 = distance(p, bh0);
        // auto R2 = distance(p, bh1);
        auto V1 = 1.0 / sqrt(2.0 * (R1 / rs - 1.0));
        // auto V2 = sqrt(params.M*params.G / (R2 - rs));

        return xzortho(b1) / length(b1) * V1;

        //auto q = R1 / (R1+R2);
        //return xzortho(b1)/length(b1)*V1 * (1.-q) + xzortho(b2)/length(b2)*V2 * q;
    };

    auto trace_block = [&](int bx, int by) {
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(-0.5, 0.5);
        std::uniform_real_distribution<> disbl(0, blocksize);

        Vec<float> tmpxyz[(50 + blocksize) * (50 + blocksize)];
        Vec<float> color;
        float color_alpha;

        const double block_div = 4.0;
        for (double bx0 = 0.5 / block_div; bx0 < blocksize; bx0 += 1.0 / block_div) {
            for (double by0 = 0.5 / block_div; by0 < blocksize; by0 += 1.0 / block_div) {
                //double x00 = disbl(gen), y00 = disbl(gen);
                double x0 = bx0;
                double y0 = by0;

                //Multisampling
                //int nr = 0;
                //while(nr < 1)
                //{
                //generate pixel of view:
                //double x0 = x00+dis(gen), y0 = y00+dis(gen);

                //Horizontal, Vertical view angle
                auto wx = (bx * blocksize + x0 - tpar.w / 2) * da * deg;
                auto wy = (by * blocksize + y0 - tpar.h / 2) * da * deg;

                //Initial position, velocity and direction vectors
                auto r0 = Vec<double>{0.0, params.ycam, params.zcam};
                auto nr0 = normalize(r0);
                auto camright = normalize(cross(up, nr0));
                auto camlook = -1.0 * nr0;
                auto camup = normalize(cross(camright, camlook));

                auto v00tmp = rot_around(camlook, camright, wy);
                auto v00 = normalize(rot_around(v00tmp, camup, -wx));
                //auto v00 = Vec<double>{sin(wx), sin(wy),     -1.0};
                //auto nv0 = normalize(v00);


                //Direction vectors for changing into orbital plane determined by the radius and velocity vectors
                auto n = normalize(cross(nr0, v00)); //normal of movement plane
                auto d = nr0;                        //radial direction
                // auto d0 = normalize(cross(up, n));   //radial direction in the plane of the disk
                auto o = normalize(cross(n, d));     //transverse direction

                //Initial values:
                auto R = length(r0);
                const auto Rlim = R * 1.1;
                auto e0 = length(v00) * sqrt(1.0 - rs / R);//energy

                auto Rdot0 = dot(nr0, v00);//radial velocity
                auto h0 = length(cross(r0, v00));//angular momentum

                auto Rdot = Rdot0;
                auto lambda0 = 0.07;//affine parameter step

                int k = 0;
                int klim = (int) (100000);//step limiter
                color = Vec<float>{0.f, 0.f, 0.f};
                color_alpha = 1.0f;
                bool hit = false;

                //calculate phi0 (initial angle, measured from the disk):
                //auto phi0 = 0.0;
                /*{
                    auto q1 = std::abs(dot(r0,  R * d));
                    auto q2 = std::abs(dot(r0, -R * d));

                    phi0 = q1 < q2 ? angle(r0, R*d) : angle(r0, -R*d);
                }*/

                auto phi = 0.0;
                auto last = r0;
                auto dl = 0.0;
                Vec<double> rrn;
                auto Rcum = 0.0;
                double Flux = 1.0;
                while (R < Rlim && k < klim) {
                    auto Rdotl = Rdot;
                    double lambda = lambda0;
                    auto rho = sqrt(1.0 - rs / R);

                    if (R < rs) {
                        k = klim;
                        hit = true;
                        break;
                    }

                    {
                        auto rr = rot_around(r0, n, phi);
                        auto nrr = normalize(rr);
                        last = rrn;
                        rrn = nrr * R;
                        auto x = rrn.x;
                        auto y = rrn.y;
                        auto z = rrn.z;
                        auto plane_sqr = sq(x) + sq(z) + sq(y);
                        dl = length(rrn - last);
                        Rcum += dl;

                        if (R < 3.0 * rs) { lambda /= 1.75; }
                        else if (R < 2.5 *
                                     rs) { lambda /= 4.0; /*if(abs(lasty) < 0.1*rs){ lambda /= 55.0; }else{ lambda /= 25.0; }*/ }

                        if (sq(rdin) <= plane_sqr && plane_sqr <= sq(rdout) && sq(y) <= sq(hdisk)) {
                            lambda = lambda0 / 9.0;

                            //sample disk density:
                            auto disk_angle = (int) ((Pi + atan2(x, z)) / (2.0 * Pi) * diskHeight);
                            auto disk_rad = (int) ((sqrt(plane_sqr) - rdin) / (rdout - rdin) * diskWidth);
                            if (disk_rad >= diskWidth) { disk_rad = diskWidth - 1; }
                            if (disk_angle >= diskHeight) { disk_angle = diskHeight - 1; }
                            disk_rad = diskWidth - 1 - disk_rad;
                            auto density = 1.0f * (float) dl * disk[disk_rad * diskHeight + disk_angle] *
                                           (1.0f - (float) abs(y / hdisk));

                            auto dvel = vel_at(Vec<double>{x, y, z});
                            auto sqmag = sqlength(dvel);
                            // auto mag = sqrt(sqmag);

                            auto gamma = 1.0 / sqrt(1.0 - sqmag);

                            auto T = /*mag * */params.T0disk;

                            auto kTtmp = normalize(rot_around(o, n, phi));
                            auto kRtmp = normalize(rot_around(d, n, phi));
                            auto KT = kTtmp * h0 / R;
                            auto KR = kRtmp * Rdot;
                            auto K = normalize(KT + KR);
                            //Boost back is missing...
                            auto view_angle_corr = sin(angle(K, Vec<double>{K.x, 0.0, K.z}));

                            auto q = dot(dvel, KT) / (gamma * rho);

                            //Doppler factor = z + 1 due to motion
                            auto D = gamma * (1 - q);
                            //Doppler factor due to gravitational well
                            auto Dz = rho_obs / rho;

                            Flux = 1.0;
                            {
                                auto rm = R / params.M;
                                auto S = sqrt(rm);
                                auto sqrt3 = sqrt(3.0);
                                auto sqrt6 = sqrt(6.0);
                                Flux = 3.0 * params.M * 1.0 / (8.0 * Pi * (rm - 3.0) * S * sq(rm)) * (S - sqrt6 +
                                                                                                      sqrt3 / 3.0 *
                                                                                                      std::log(
                                                                                                          (S + sqrt3) *
                                                                                                          (sqrt6 -
                                                                                                           sqrt3) /
                                                                                                          ((S - sqrt3) *
                                                                                                           (sqrt6 +
                                                                                                            sqrt3))));
                            }

                            color = color +
                                    color_alpha * black_body_xyz((float) (T * (Dz * D))) * Flux * view_angle_corr *
                                    density * 2.0e1f / sq((float) Rcum);

                            color_alpha *= /*(1.0f - 4.0f*dl/hdisk)**/0.99965f * (1.0f - sqrt(density));//0.99

                            //color = color + color_alpha*black_body_xyz( (float)(T) )*Flux;//*2.0e1f/sq((float)Rcum)*view_angle_corr;

                            if (color_alpha <= 0.0001f) {
                                color_alpha = 0.0001f;
                                hit = true;
                                break;
                            }
                        }

                    }

                    phi += lambda * h0 / sq(R);
                    Rdot += lambda * (params.M / sq(R * rho) * (sq(Rdotl) - sq(e0)) + sq(h0 * rho) / cube(R));
                    R += lambda * Rdotl;

                    k += 1;
                }
                if (!hit) {
                    auto rr = rot_around(r0, n, phi);
                    auto nrr = normalize(rr);
                    auto rrn = nrr * R;
                    auto x = rrn.x;
                    auto y = rrn.y;
                    auto z = rrn.z;

                    /*auto x = P.x;
                    auto y = P.y;
                    auto z = P.z;*/

                    auto r = sqrt(sq(x) + sq(y) + sq(z));
                    x /= r;
                    y /= r;
                    z /= r;

                    // auto theta_tmp = acos(z);
                    // auto phi_tmp = Pi + atan2(y, x);
                    auto tr = (int) (acos(z) / deg);
                    auto pr = 180 + (int) (atan2(y, x) / deg);
                    auto mi = pr * 180 + tr;
                    double dmin = 100000.0;
                    int ss = -1;
                    auto const& starblock = starmap[mi];
                    for (int s = 0; s < (int) starblock.size(); ++s) {
                        auto const& star = stars[starblock[s]];
                        auto theta = (double) star.angles.theta;
                        auto phi = (double) star.angles.phi;
                        auto xs = sin(theta) * cos(Pi + phi);
                        auto ys = sin(theta) * sin(Pi + phi);
                        auto zs = cos(theta);

                        auto d = sq(x - xs) + sq(y - ys) + sq(z - zs);
                        if (d < dmin) {
                            ss = s;
                            dmin = d;
                        }
                    }

                    if (ss >= 0) {
                        auto const& star = stars[starblock[ss]];
                        if (dmin < star.d) {
                            if (color_alpha > 0.0001f) {
                                color = color + 1e-12f * pow(color_alpha, 14.0f) * star.colxyz;
                            } else {
                                color = star.colxyz;
                            }
                        }
                    }
                }

                auto Lo = [](auto const& x, auto const& w) { return 1.0 / (w * (1.0 + sq(x / w))); };
                auto Bump = [](auto const& x, auto const& w) {
                    if (abs(x) > w) { return 0.0; }
                    return exp(-1.0 / (1.0 + sq(abs(x) / w)));
                };

                if (length(color) != 0.0f) {
                    // auto cs = (color.x+color.y+color.z);
                    static const int L = 24;
                    for (int yi = -L; yi <= L; ++yi) {
                        for (int xi = -L; xi <= L; ++xi) {
                            auto xp = (int) (25 + x0 + xi);
                            auto yp = (int) (25 + y0 + yi);
                            auto Ex = (float) (sq(xi) + sq(yi));
                            //auto A = (1.0 + pow(-Ex*(cs / 5e-3), 0.2));
                            auto dec = exp(-0.5 * Ex / sq(24.0 / 5.5)) * Bump(sqrt(Ex),
                                                                              24.0f);//(1.0f - tanh(Ex-10.0f))*(1.0f - sq(sqrt(Ex)/20.0f));//exp(Ex / sq(4.0));
                            auto ratx = (float) (dec * Lo(-(Ex) / 2.0f, 0.1 * 0.74));//(A*exp(Ex/sq(0.25*0.74)));
                            auto raty = (float) (dec * Lo(-(Ex) / 2.0f, 0.1 * 0.77));//(A*exp(Ex/sq(0.25*0.77)));
                            auto ratz = (float) (dec * Lo(-(Ex) / 2.0f, 0.1 * 0.89));//(A*exp(Ex/sq(0.25*0.89)));
                            if (xp >= 0 && yp >= 0 && xp < blocksize + 50 && yp < blocksize + 50) {
                                auto c0 = tmpxyz[yp * (blocksize + 50) + xp];
                                c0.x += color.x * ratx;
                                c0.y += color.y * raty;
                                c0.z += color.z * ratz;
                                tmpxyz[yp * (blocksize + 50) + xp] = c0;
                            }
                        }
                    }

                }
                //nr += 1;
                //} //Multisample
            }//block y
        }//block x

        int y0 = by * blocksize - 25;
        int x0 = bx * blocksize - 25;
        {
            std::lock_guard<std::mutex> lock(mxyz);
            for (int yi = 0; yi < blocksize + 50; ++yi) {
                for (int xi = 0; xi < blocksize + 50; ++xi) {
                    auto xp = x0 + xi;
                    auto yp = y0 + yi;
                    if (xp >= 0 && yp >= 0 && xp < tpar.w && yp < tpar.h) {
                        auto c0 = xyz[yp * tpar.w + xp];
                        c0 = c0 + tmpxyz[yi * (blocksize + 50) + xi];
                        xyz[yp * tpar.w + xp] = c0;
                    }
                }
            }
        }
    };//trace block

    std::vector<std::future<void>> threads(nths);
    for (auto& th : threads) {
        th = std::async(std::launch::async, trace_block, tpar.bx, tpar.by);
        tpar.step();
    }
    for (auto& th : threads) { (void) th.get(); }
}

std::mt19937 gen(42);

auto gen_random_value_in_interval = [](auto const& lo, auto const& hi) {
    std::uniform_real_distribution<> dis(lo, hi);
    return dis(gen);
};

auto gen_random_value_around = [](auto const& mu, auto const& sigma) {
    std::lognormal_distribution<> dis(mu, sigma);
    return dis(gen);
};

auto gen_random_value_exp = [](auto const& lambda) {
    std::exponential_distribution<> dis(lambda);
    return dis(gen);
};

auto gen_random_theta_phi = [] {
    std::uniform_real_distribution<> dis(0.0, 1.0);
    return ThetaPhi<double>{std::acos(1.0 - 2.0 * dis(gen)), 2.0 * Pi * dis(gen)};
};

#if _WIN32
int wWinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, PWSTR pCmdLine, int nCmdShow)
#else

int main(int argc, char* argv[])
#endif
{
    int width = 1920, height = 1080;
    sf::RenderWindow window(sf::VideoMode(width, height), "Gravitation Raytrace");

    sf::RectangleShape rect(sf::Vector2f(width * 1.0f, height * 1.0f));
    rect.setFillColor(sf::Color(255, 255, 255));
    rect.setOutlineColor(sf::Color(255, 255, 255));

    sf::Font font;
#if _WIN32
    if (!font.loadFromFile("C:\\Windows\\Fonts\\arial.ttf")){ return -1; }
#else
    if (!font.loadFromFile("/usr/share/fonts/truetype/msttcorefonts/arial.ttf")) { return -1; }
#endif

    sf::Text text;
    text.setFont(font); // font is a sf::Font
    text.setPosition(30, 30);
    text.setString("0");
    text.setCharacterSize(12);
#if _WIN32
    text.setFillColor(sf::Color(0, 255, 0));
#else
    text.setColor(sf::Color(0, 255, 0));
#endif


    std::vector<Star> stars(200000);
    std::vector<std::vector<size_t>> starmap(360 * 180);

    for (size_t i = 0; i < stars.size(); ++i) {
        float T;
        Vec<float> bb;
        auto dg = 1.0;
        bool accepted = false;
        while (!accepted) {
            const auto Tmag = 3000.0f;
            const auto spread = 0.7f;
            auto TI = 1.0f;
            auto giant = gen_random_value_in_interval(0.0, 100.0);
            auto var = 20000.0f * pow(10.0f, (float) gen_random_value_in_interval(0.5f, 10.5f));
            dg = 1.0;
            if (giant > 98.9) {
                //Giant
                T = (float) (Tmag * gen_random_value_around(0.0, spread));
                TI = 1e2f * var;
                dg = 2.2;
            } else if (giant > 95.0) {
                T = (float) (Tmag * gen_random_value_around(0.0, spread));
                TI = 5e1f * var;
                dg = 1.8;
            } else if (giant > 92.0) {
                T = (float) (Tmag * gen_random_value_around(0.0, spread));
                TI = 1e1f * var;
                dg = 1.25;
            } else {
                //regular star
                T = (float) (Tmag * gen_random_value_around(0.0, spread));
                TI = var;
            }

            auto x = gen_random_value_in_interval(0.1f, 1500000.f);
            auto y = gen_random_value_in_interval(0.1f, 1500000.f);
            auto z = gen_random_value_in_interval(0.1f, 1500000.f);
            auto sqr = sq(x) + sq(y) + sq(z);

            bb = black_body_xyz(T) * TI;
            if (gen_random_value_in_interval(-1.0, 1.0) > 0.0) {
                //binary partner
                auto T2 = (float) (Tmag * gen_random_value_around(0.0, spread));
                bb = bb + black_body_xyz(T2) * pow(10.0f, (float) gen_random_value_in_interval(0.5f, 10.5f));
            }
            bb = bb / sqr;
            if (bb.sum() > 1e-8) { accepted = true; }
        }

        auto as = gen_random_theta_phi();
        stars[i].angles.theta = (float) as.theta;
        stars[i].angles.phi = (float) as.phi;

        auto th = (int) (as.theta / deg);
        auto ph = (int) (as.phi / deg);

        auto helper = [&](int t, int p, size_t idx) mutable {
            if (p >= 0 && p <= 359 && t >= 0 && t <= 179) {
                starmap[p * 180 + t].push_back(idx);
            }
        };
        helper(th, ph - 1, i);
        helper(th, ph, i);
        helper(th, ph + 1, i);
        helper(th + 1, ph - 1, i);
        helper(th + 1, ph, i);
        helper(th + 1, ph + 1, i);
        helper(th - 1, ph - 1, i);
        helper(th - 1, ph, i);
        helper(th - 1, ph + 1, i);

        //arcdiameter
        stars[i].d = (float) sq(dg * 2.5e-4 * gen_random_value_around(0.0, 0.25));
        stars[i].T = T;
        stars[i].colxyz = bb;
    }

    int diskW = diskHeight;
    int diskH = diskWidth;
    std::vector<float> disk(diskW * diskH, 0.0f);

    for (int l = 2; l <= 8; ++l) {
        auto f = pow(2.0, 1.0 * l);
        // auto rf = 1.0/f;
        for (int n = 0; n < 125 * f; ++n) {
            auto y0 = (int) (sq(sq(gen_random_value_in_interval(0.0, 1.0))) * (diskH - 15));
            auto h = (int) gen_random_value_in_interval(1, 2 + f / 150.);
            auto x0 = gen_random_value_in_interval(0, diskW);
            auto w = gen_random_value_in_interval(/*0.02*diskW*/0.0, 0.1 * diskW * (0.1 + (double) y0 / (1.0 * diskH)));
            for (int yn = y0; yn < y0 + h; ++yn) {

                auto x = (int) x0;
                for (int xn = 0; xn < w; ++xn) {
                    if (x == diskW) { x = 0; }
                    disk[yn * diskW + x] += (float) (sq(f) * sq(y0 / (1.0 *
                                                                      diskH)));//(float)(f*y0*(diskH-y0)/(1.0f*sq(diskH)));
                    x += 1;
                }
            }
        }
    }

    auto diskImax = *std::max_element(disk.cbegin(), disk.cend());
    for (auto& x:disk) { x = cube(x / diskImax); }

    Params params;
    params.G = 1.0;
    params.M = 1.0;
    params.fovw = 70.0;
    params.zcam = 32.0;
    params.ycam = -1.95;//-1.95;
    params.T0disk = 9000.0;

    float Ifactor = 1.0;//133.657e-6;
    float Iexp = 0.5f;//0.5f;
    float minI, maxI;

    bool trace_stop = false;
    bool image_saved = false;
    int image_id = 0;

    TraceParams tpar;
    tpar.finish = &trace_stop;

    std::vector<Vec<float>> xyz_image;
    std::vector<sf::Color> result_image;

    sf::Image img;
    sf::Texture tex;
    sf::Sprite sprite;

    auto resize_img = [&] {
        tpar.w = width;
        tpar.h = height;
        tpar.reset();
        xyz_image.resize(width * height);
        result_image.resize(width * height);
        for (int y = 0; y < height; y += 1) {
            for (int x = 0; x < width; x += 1) {
                xyz_image[y * width + x] = Vec<float>{0.f, 0.f, 0.f};
                result_image[y * width + x] = sf::Color::Black;
            }
        }
    };
    resize_img();

    auto update_img = [&] {
        if (!trace_stop) {
            trace_n(tpar, stars, starmap, xyz_image, result_image, disk, params);
        }

        minI = 1e19f;
        maxI = std::accumulate(xyz_image.cbegin(), xyz_image.cend(), 0.0f, [&](auto const& acc, auto const& col) {
            auto I = col.y;
            minI = I > 0 ? std::min(minI, I) : minI;
            return std::max(acc, I);
        });
        maxI *= Ifactor;
        //minI *= 1.e1f;

        auto lmaxI = maxI;
        auto lminI = minI;
        auto range = lmaxI - lminI;
        for (int y = 0; y < height; y += 1) {
            for (int x = 0; x < width; x += 1) {
                auto xyz = xyz_image[y * width + x];
                auto q = (float) pow(xyz.y / range, Iexp);
                if (q < 0.0f) { q = 0.0f; };
                auto rgb = convert_xyz_to_rgb(xyz);
                rgb = gamma_correct(redistribute_rgb(rgb / (1.0f) * q), 2.4f);
                result_image[y * width + x] = rgb.sfColor();
            }
        }

        img.create(width, height, (sf::Uint8*) result_image.data());
        tex.loadFromImage(img);
        sprite.setTexture(tex, true);
        sprite.setPosition(sf::Vector2f(0, 0));

        if (trace_stop && !image_saved) {
            img.saveToFile(std::string("result") + std::to_string(image_id) + ".png");
            image_id += 1;
            image_saved = true;
        }
    };

    auto draw_text = [&] {
        if (!trace_stop) {
            text.setString(std::to_string(tpar.ratio_complete() * 100.0).substr(0, 6) + std::string("% complete"));
            text.setPosition(30, 30);
            window.draw(text);
        }
        /*text.setString(std::string("Mass = ") + std::to_string(params.M));
        text.setPosition(30, 30);
        window.draw(text);

        text.setString(std::string("z = ") + std::to_string(params.zcam));
        text.setPosition(30, 50);
        window.draw(text);

        text.setString(std::string("y = ") + std::to_string(params.ycam));
        text.setPosition(30, 70);
        window.draw(text);

        text.setString(std::string("Imax = ") + std::to_string(maxI));
        text.setPosition(30, 90);
        window.draw(text);

        text.setString(std::string("Imin = ") + std::to_string(minI));
        text.setPosition(30, 110);
        window.draw(text);

        text.setString(std::string("Ifactor = ") + std::to_string(1e6*Ifactor));
        text.setPosition(30, 130);
        window.draw(text);

        text.setString(std::string("Iexp = ") + std::to_string(Iexp));
        text.setPosition(30, 150);
        window.draw(text);*/
    };

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) { window.close(); }
                //else if(event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::Key::Left) { params.zcam -= 0.1; resize_img(); }
                //else if(event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::Key::Right){ params.zcam += 0.1; resize_img(); }
            else if (event.type == sf::Event::KeyPressed &&
                     event.key.code == sf::Keyboard::Key::Up) { Ifactor *= 1.5;/*params.ycam += 0.1; resize_img();*/ }
            else if (event.type == sf::Event::KeyPressed &&
                     event.key.code == sf::Keyboard::Key::Down) { Ifactor /= 1.5;/*params.ycam -= 0.1; resize_img();*/ }
            else if (event.type == sf::Event::KeyReleased &&
                     event.key.code == sf::Keyboard::Key::Tab) { trace_stop = !trace_stop; }
            else if (event.type == sf::Event::KeyPressed && event.key.code ==
                                                            sf::Keyboard::Key::Subtract) { Iexp -= 0.005; /*params.M -= 0.05; resize_img();*/ }
            else if (event.type == sf::Event::KeyPressed &&
                     event.key.code == sf::Keyboard::Key::Add) { Iexp += 0.005; /*params.M += 0.05; resize_img();*/ }
            else if (event.type == sf::Event::Resized) {
                width = event.size.width;
                height = event.size.height;
                window.setView(sf::View(sf::FloatRect(0.0f, 0.0f, width * 1.0f, height * 1.0f)));
                glViewport(0, 0, width, height);
                rect.setSize(sf::Vector2f(width * 1.0f, height * 1.0f));
                resize_img();
            } else if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::Key::Space) {
                img.saveToFile(std::string("result") + std::to_string(image_id) + ".png");
                image_id += 1;
            }
        }

        update_img();

        window.clear();
        window.draw(rect);
        window.draw(sprite);
        draw_text();
        window.display();
    }

    return 0;
}