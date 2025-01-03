#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <vector>

#include "fenster.h"

#define FPS 120

#define BLACK 0x000000
#define WHITE 0xFFFFFF
#define RED 0xFF0000

typedef uint8_t u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;
typedef int32_t i32;
typedef int64_t i64;
typedef float f32;
typedef double f64;

// Note: Y-axis is positive in downward direction.
typedef struct point {
    i32 x = 0;
    i32 y = 0;
} point;

// TODO AABBs on a torus are tricky to get right. Fix this.
typedef struct aabb {
    // Guarantees:
    //   0 <= xmin <= xmax
    //   0 <= ymin <= ymax
    //   xmin < width
    //   ymin < height
    // It's _not_ guaranteed that xmax < width, or ymax < height (our world is a torus).
    i32 xmin = 0;
    i32 xmax = 0;
    i32 ymin = 0;
    i32 ymax = 0;
} aabb;

typedef struct cluster {
    std::vector<point> points = {};
    // The bounding box is not guaranteed minimal. E.g., it does _not_ take into account the fact that the world is a
    // torus. E.g., if the cluster straddles both the X and Y axis, then the bounding box may be the entire world.
    aabb bounds = {};
    u32 color = WHITE;
} cluster;

typedef struct world {
    i32 w = 0;
    i32 h = 0;
    std::vector<cluster> clusters = {};
} world;

void expand_aabb(aabb *rect, point p) {
    rect->xmin = std::min(rect->xmin, p.x);
    rect->xmax = std::max(rect->xmax, p.x);
    rect->ymin = std::min(rect->ymin, p.y);
    rect->ymax = std::max(rect->ymax, p.y);
}

// TODO This makes overly-large AABBs because it doesn't consider clusters that straddle the axes.
aabb make_aabb(std::vector<point> points) {
    aabb result;
    if (points.size() == 0) {
        return result;
    }
    result = {points[0].x, points[0].x,
              points[0].y, points[0].y};
    for (int i = 1; i < points.size(); ++i) {
        expand_aabb(&result, points[i]);
    }
    return result;
}

aabb merge_aabb(const aabb a, const aabb b) {
    aabb result;
    result.xmin = std::min(a.xmin, b.xmin);
    result.xmax = std::max(a.xmax, b.xmax);
    result.ymin = std::min(a.ymin, b.ymin);
    result.ymax = std::max(a.ymax, b.ymax);
    return result;
}

// Return true if the two AABBs might be touching.
// (Note this may sometimes return true when they are not touching.)
bool touching(const aabb a, const aabb b) {
    if (a.xmin > b.xmax + 1) return false;
    if (a.ymin > b.ymax + 1) return false;
    if (b.xmin > a.xmax + 1) return false;
    if (b.ymin > a.ymax + 1) return false;
    return true;
}

// Will find the cluster at (x, y), remove it from grid, and load it into c.
// Precondition: grid[y*w + x] == 1.
cluster extract_cluster(u8 *grid, i32 w, i32 h, i32 x, i32 y) {
    cluster c;
    c.color = rand();
    if (grid[y * w + x] == 0) {
        fprintf(stderr, "Warning: extract_cluster(): The point (%d, %d) is empty.\n", x, y);
        // Nothing to do -- return empty cluster.
        return c;
    }
    // Bucket fill, on grid as torus (i.e., wraps around).
    std::vector<point> stack = {{x, y},};
    while (!stack.empty()) {
        point p = stack.back();
        stack.pop_back();
        if (grid[p.y * w + p.x] != 0) {
            grid[p.y * w + p.x] = 0;
            c.points.push_back(p);
            point adj[4];
            adj[0] = {(p.x+w-1)%w, p.y};
            adj[1] = {(p.x+w+1)%w, p.y};
            adj[2] = {p.x, (p.y+h+1)%h};
            adj[3] = {p.x, (p.y+h-1)%h};
            for (auto &a : adj) {
                stack.push_back(a);
            }
        }
    }
    // Superfluous for the time being.
    //c.bounds = make_aabb(c.points);
    return c;
}

void destroy(world *wld) {
    wld->w = 0;
    wld->h = 0;
    wld->clusters.clear();
}

void populate(world *wld, i32 w, i32 h, f64 density) {
    wld->w = w;
    wld->h = h;

    u8 *grid;
    grid = (u8 *) malloc(sizeof(u8) * w * h);

    // Randomly create matter.
    for (int i = 0; i < h; ++i) {
        for (int j = 0; j < w; ++j) {
            // TODO Improve RNG
            //   - our own RNG, e.g., JSF
            //   - option to populate combinatorially (n choose k), with density exactly as given
            grid[i * w + j] = rand() <= (f64)RAND_MAX * density;
        }
    }

    // Find initial clusters.
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            if (grid[y * w + x]) {
                cluster c = extract_cluster(grid, w, h, x, y);
                wld->clusters.push_back(c);
            }
        }
    }

    free(grid);
    grid = NULL;
    return;
}

u32 largest_cluster_size(world *wld) {
    u32 largest_size = 0;
    for (const auto &c : wld->clusters) {
        largest_size = std::max(largest_size, (u32)c.points.size());
    }
    return largest_size;
}

u32 num_clusters(world *wld) {
    return (u32)wld->clusters.size();
}

u32 total_mass(world *wld) {
    size_t result = 0;
    for (const auto &c : wld->clusters) {
        result += c.points.size();
    }
    return (u32)result;
}

// L1 distance: Sum of X and Y distances, taking into account grid wrap-around (e.g., the distance between NW corner and
// SE corner of grid is 2).
inline i32 distance(const point &a, const point &b, i32 w, i32 h) {
    i32 dist_x = abs(b.x - a.x);
    dist_x = std::min(dist_x, w - dist_x);
    i32 dist_y = abs(b.y - a.y);
    dist_y = std::min(dist_y, h - dist_y);
    return dist_x + dist_y;
}

bool contact(const cluster &a, const cluster &b, i32 w, i32 h) {
    if (!touching(a.bounds, b.bounds)) {
        return false;
    }
    for (const auto &p1 : a.points) {
        for (const auto &p2 : b.points) {
            if (distance(p1, p2, w, h) <= 1) {
                return true;
            }
        }
    }
    return false;
}

void evolve(world *wld) {
    std::vector<cluster> &clusters = wld->clusters;
    for (auto &c : wld->clusters) {
        // Need this instead of 0 to deal with C modulo operator shenanigans.
        int dx = wld->w;
        int dy = wld->h;
        switch (rand() % 4) {
        default:
        case 0: dx += 1; break;
        case 1: dx -= 1; break;
        case 2: dy += 1; break;
        case 3: dy -= 1; break;
        }
        for (auto &p : c.points) {
            p.x = (p.x + dx) % wld->w;
            p.y = (p.y + dy) % wld->h;
        }
        // TODO Buggy: When AABBs move across X/Y axis, this can break.
        //c.bounds.xmin += dx;
        //c.bounds.xmax += dx;
        //c.bounds.ymin += dy;
        //c.bounds.ymax += dy;
        // For now, just re-generate AABBs at every step -- this is quite slow, but better than nothing.
        c.bounds = make_aabb(c.points);
    }

    // Check clusters for collisions, and merge if needed.
    u32 nclusters = clusters.size();
    for (u32 j = 0; j < nclusters; ++j) {
        for (u32 k = j+1; k < nclusters; ++k) {
            cluster &cj = clusters[j];
            cluster &ck = clusters[k];
            if (contact(cj, ck, wld->w, wld->h)) {
                if (cj.points.size() < ck.points.size()) {
                    // Keep the color of the larger cluster.
                    // TODO Bug: This may result in the wrong color if more than two clusters combine in a single
                    // evolve() step.
                    cj.color = ck.color;
                }
                cj.points.insert(cj.points.begin(), ck.points.begin(), ck.points.end());
                // TODO AABB merging.
                //cj.bounds = merge_aabb(cj.bounds, ck.bounds);
                clusters.erase(clusters.begin() + k);
                --nclusters;
            }
        }
    }
    // TODO Come up with a better way to deal with overlapping (clusters running into and over each other) -- it
    // should not destroy mass or have "overlapping" mass.
}

void render(const world *wld, const fenster *f) {
    // Clear display.
    for (int i = 0; i < f->width; i++) {
        for (int j = 0; j < f->height; j++) {
            fenster_pixel(f, i, j) = BLACK;
        }
    }
    // Paint clusters.
    for (const auto &c : wld->clusters) {
        u32 color = c.color;
        color |= 0x404040;  // For visibility.
        for (const auto &p : c.points) {
            fenster_pixel(f, p.x, p.y) = color;
        }
    }
}

i32 run(world *wld, f64 threshold, bool display = true, bool verbose = false) {
    u32 *buf = (u32 *)malloc(sizeof(u32) * wld->w * wld->h);
    struct fenster f = {
        .title = "Aggregation",
        .width = wld->w,
        .height = wld->h,
        .buf = buf,
    };
    if (display) {
        fenster_open(&f);
    }

    u64 step = 0;
    u32 agglom = 0;
    u32 agglom_target = total_mass(wld) * threshold;
    bool completed = false;
 
    i64 prev_render = 0;
    if (display) {
        fenster_sleep(1000/FPS); // Fix display bug: no display for first frame.
        prev_render = fenster_time();
        render(wld, &f);
    }
    while (true) {
        if (display && fenster_loop(&f) != 0) {
            break;
        }
        agglom = largest_cluster_size(wld);
        if (verbose) {
            fprintf(stdout,
                    "Step %lu: Clusters: %u. Largest cluster size: %u.\n",
                    step, num_clusters(wld), agglom);
        }
        if (agglom >= agglom_target) {
            completed = true;
            break;
        }

        evolve(wld);
        ++step;
        if (display) {
            i64 now = fenster_time();
            if (now - prev_render > 1000/FPS) {
                prev_render = now;
                render(wld, &f);
            }
            //fenster_sleep(500);
        }
    }

    if (display) {
        fenster_close(&f);
    }
    if (completed) {
        fprintf(stdout,
                "Achieved %.2f%% agglomeration in %lu steps.\n",
                100 * (f64)agglom / (f64)(total_mass(wld)),
                step);
    }
    else {
        fprintf(stdout, "Aborted after %lu steps.\n", step);
    }

    free(buf);
    buf = NULL;
    return step;
}

int main(int argc, char **argv) {
    if (argc != 6) {
        fprintf(stderr, "Usage: aggregation <width> <height> <density> <threshold> <trials>\n");
        fprintf(stderr, "  width, height: Size of grid, e.g., 100 by 100\n");
        fprintf(stderr, "  density: Initial proportion of cells with matter, e.g., 0.1. Must be between 0 and 1.\n");
        fprintf(stderr, "  threshold: The target mass proportion for the largest cluster, e.g., 0.9. Must be between 0 and 1.\n");
        fprintf(stderr, "  trials: The number of simulation runs to perform.\n");
        return 1;
    }
    i32 w = atoi(argv[1]);
    i32 h = atoi(argv[2]);
    f64 density = atof(argv[3]);
    f64 threshold = atof(argv[4]);
    i32 trials = atoi(argv[5]);
    bool display = true;
    if (trials > 1) {
        fprintf(stdout, "Running multiple trials -- hiding output.\n");
        display = false;
    }

    if (trials <= 0 || w <= 0 || h <= 0) {
        fprintf(stderr, "Invalid parameters.\n");
    }
    if (density < 0.0 || density >= 1.0) {
        fprintf(stderr, "Invalid density.\n");
    }
    if (threshold < 0.0 || threshold >= 1.0) {
        fprintf(stderr, "Invalid threshold.\n");
    }


    srand(time(NULL));

    world wld;
    i32 min_steps = 0;
    i32 max_steps = 0;
    f64 mean_steps = 0.0;
    f64 sum_sqr_steps = 0.0;
    f64 stddev_steps = 0.0;
    for (i32 i = 0; i < trials; ++i) {
        printf("Running trial %d/%d...\n", i, trials);
        populate(&wld, w, h, density);
        i32 steps = run(&wld, threshold, display, false);
        destroy(&wld);
        min_steps = (i == 0) ? steps : std::min(min_steps, steps);
        max_steps = std::max(max_steps, steps);
        mean_steps += (f64)steps / trials;
        sum_sqr_steps += steps*steps;
    }
    stddev_steps = sqrt((f64)sum_sqr_steps/trials - mean_steps*mean_steps);
    printf(
        "Ran %d trials (size %d×%d, density %.3f, threshold %.3f). Steps (min, max, mean±stddev): %d, %d, %.1f ± %.0f.\n",
        trials, w, h, density, threshold,
        min_steps, max_steps, mean_steps, stddev_steps);

    return 0;
}
