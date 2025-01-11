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

// TODO Unit tests, and end-to-end tests.

// Note: Y-axis is positive in downward direction.
typedef struct point {
    i32 x = 0;
    i32 y = 0;
} point;

// A "Toroidal AABB" (Axis-Aligned Bounding Box) comprises a "toroidal interval" for each axis.
// The semantics for a toroidal interval is as follows, where N is width (for x axis) or height (for y axis).
//   It's always the case that min and max satisfy 0 <= min < N and 0 <= max < N.
//   if (min == max): All points x in body satisfy x == min (== max).
//   if (min < max): All points x in body satisfy min <= x <= max.
//            (This includes the case 0 == min && max == width - 1, which is the trivial bounding interval.)
//   if (max < min): All points x satisfy either 0 <= x <= max or min <= x < width.
//            (This is the "wraparound" case, where the body straddles 0 for this axis.)
typedef struct taabb {
    i32 xmin = 0;
    i32 xmax = 0;
    i32 ymin = 0;
    i32 ymax = 0;
} taabb;

typedef struct cluster {
    std::vector<point> points = {};
    // The bounding box is not guaranteed minimal. E.g., it does _not_ take into account the fact that the world is a
    // torus. E.g., if the cluster straddles both the X and Y axis, then the bounding box may be the entire world.
    taabb bounds = {};
    u32 color = WHITE;
} cluster;

typedef struct world {
    i32 w = 0;
    i32 h = 0;
    std::vector<cluster> clusters = {};
} world;


// Given a bitmap, set (min, max) to a toroidal interval that contains all 1's in the bitmap.  The interval is
// guaranteed to contain all 1's, and is guaranteed to be minimal if the bitmap's 1's are contiguous (i.e., connected in
// the 1-torus).
void make_taabb_1D(const std::vector<u8> &bitmap, i32 *min, i32 *max) {
    bool started = false;
    bool ended = false;
    u32 n = bitmap.size();
    *min = 0;
    *max = 0;
    for (int x = 0; x < n; ++x) {
        if (!started) {
            if (bitmap[x]) {
                *min = x;
                *max = n - 1;
                started = true;
            }
        }
        else {  // started
            if (!ended) {
                if (!bitmap[x]) {
                    *max = x - 1;
                    ended = true;
                }
            }
            else {  // ended
                if (bitmap[x]) {
                    // Wraparound for this axis.
                    *min = x;
                    break;
                }
            }
        }
    }
}

// Return a TAABB containing all given points. The TAABB is guaranteed to be minimal, provided that (a) points is
// nonempty, and (b) points is contiguous (i.e., connected in the torus).
taabb make_taabb(std::vector<point> points, i32 w, i32 h) {
    if (points.size() == 0) {
        // We don't allow empty TAABBs, so we return the smallest-possible nonempty TAABB.
        return { 0, 0, 0, 0 };
    }
    if (points.size() == 1) {
        // Save time.
        return { .xmin = points[0].x, .xmax = points[0].x,
                 .ymin = points[0].y, .ymax = points[0].y };
    }

    // Bitmaps.
    std::vector<u8> xs(w, 0);
    std::vector<u8> ys(h, 0);
    for (const point &p : points) {
        xs[p.x] = 1;
        ys[p.y] = 1;
    }

    taabb result;
    make_taabb_1D(xs, &result.xmin, &result.xmax);
    make_taabb_1D(ys, &result.ymin, &result.ymax);
    return result;
}

// Test whether two toroidal intervals touch.
// See definition of TAABB for specification of toroidal interval.
// Parameters:
//   a_min, a_max: First toroidal interval.
//   b_min, b_max: Second toroidal interval.
//   n: Size of torus, [0, ..., n). Must be greater than or equal to 1.
// Return:
//   true: the intervals may be touching (i.e., there may be two points of distance 0 or 1).
//   false: the intervals are definitely not touching (i.e., distance is at least 2).
bool touching_1D(i32 a_min, i32 a_max, i32 b_min, i32 b_max, i32 n) {
    if (a_min <= a_max && b_min <= b_max) {
        // Neither interval wraps around.
        if (b_max + 1 < a_min) return (b_min == 0) && (a_max == n-1);
        if (a_max + 1 < b_min) return (a_min == 0) && (b_max == n-1);
        return true;
    }
    else if (a_min > a_max && b_min <= b_max) {
        // Only interval A wraps around.
        return (b_min <= a_max + 1) || (a_min - 1 <= b_max);
    }
    else if (a_min <= a_max && b_min > b_max) {
        // Only interval B wraps around.
        return (a_min <= b_max + 1) || (b_min - 1 <= a_max);
    }
    else {
        // Both intervals wrap around.
        return true;
    }
}

// Return true if the two TAABBs might be touching.
// (This may sometimes return true even when they are not touching.)
bool touching(const taabb a, const taabb b, i32 w, i32 h) {
    return
        touching_1D(a.xmin, a.xmax, b.xmin, b.xmax, w) &&
        touching_1D(a.ymin, a.ymax, b.ymin, b.ymax, h);
}

// Return true if the point might be touching (i.e., adjacent to or overlapping) the TAABB.
// (Note this may sometimes return true when they are not touching.)
bool touching(const point p, const taabb rect, i32 w, i32 h) {
    return
        touching_1D(p.x, p.x, rect.xmin, rect.xmax, w) &&
        touching_1D(p.y, p.y, rect.ymin, rect.ymax, h);
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
    c.bounds = make_taabb(c.points, w, h);
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
            //   - option to provide seed on command line
            //   - option to populate combinatorially (n choose k), with density _exactly_ as given.
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
    // Do some preliminary checks to save work (otherwise, the time is quadratic in cluster size).

    if (!touching(a.bounds, b.bounds, w, h)) {
        return false;
    }

    std::vector<point> a_test_points;
    std::vector<point> b_test_points;
    for (const auto &p : a.points) {
        if (touching(p, b.bounds, w, h)) {
            a_test_points.push_back(p);
        }
    }
    for (const auto &p : b.points) {
        if (touching(p, a.bounds, w, h)) {
            b_test_points.push_back(p);
        }
    }

    for (const auto &p1 : a_test_points) {
        for (const auto &p2 : b_test_points) {
            if (distance(p1, p2, w, h) <= 1) {
                return true;
            }
        }
    }

    return false;
}

void evolve(world *wld) {
    std::vector<cluster> &clusters = wld->clusters;
    u32 nclusters = clusters.size();

    // Move each cluster, checking for collisions and merging.
    //u32 surplus_motion = 0;
    for (u32 j = 0; j < nclusters; ++j) {
        cluster *cj = &clusters[j];

        // We need this instead of 0 to deal with C modulo operator shenanigans.
        int dx = wld->w;
        int dy = wld->h;
        switch (rand() % 4) {
        default:
        case 0: dx += 1; break;
        case 1: dx -= 1; break;
        case 2: dy += 1; break;
        case 3: dy -= 1; break;
        }
        for (auto &p : cj->points) {
            p.x = (p.x + dx) % wld->w;
            p.y = (p.y + dy) % wld->h;
        }

        // Slide the TAABB (this is faster than calling make_taabb()).
        cj->bounds.xmin += dx;
        cj->bounds.xmax += dx;
        cj->bounds.ymin += dy;
        cj->bounds.ymax += dy;
        cj->bounds.xmin %= wld->w;
        cj->bounds.xmax %= wld->w;
        cj->bounds.ymin %= wld->h;
        cj->bounds.ymax %= wld->h;

        // Check for collisions with other clusters, and merge if needed.
        // TODO Speed this up: Use a spacial data structure so we don't have to check n^2 pairs for collision.

        // Note: It's simpler to check for collisions here, rather than first moving all the clusters and then checking
        // for collisions. If we checked for collisions after moving all clusters, we could end up with "overlap" where
        // two clusters move onto one another simultaneously, in which case we'd need to undo one or more translations.
        for (u32 k = 0; k < nclusters; ++k) {
            // Clusters of index k < j have already moved.
            // Cluster of index k == j is the one we've moved just now.
            // Clusters of index k > j have yet to be moved.
            if (k == j) {
                continue;
            }
            cluster *ck = &clusters[k];
            if (contact(*cj, *ck, wld->w, wld->h)) {
                if (cj->points.size() < ck->points.size()) {
                    // Keep the color of the larger cluster.
                    cj->color = ck->color;
                }
                cj->points.insert(cj->points.end(), ck->points.begin(), ck->points.end());

                // TODO Smarter TAABB merging.
                //cj->bounds = merge_taabb(cj->bounds, ck->bounds);
                cj->bounds = make_taabb(cj->points, wld->w, wld->h);

                // Note that this wastes cluster k's thermal energy for this step when k > j.
                clusters.erase(clusters.begin() + k);
                if (k < j) {
                    // Index moved; repair it.
                    --j;
                    cj = &clusters[j];
                }
                //else {
                //    // Cluster k has yet to move, but we've already merged it into cluster j. "Save" the energy,
                //    // and move cluster j again.
                //    ++surplus_motion;
                //}
                --k;
                --nclusters;
            }
        }
        //if (surplus_motion > 0) {
        //    --surplus_motion;
        //    --j;
        //}
    }

    // The old way: Move all clusters first, and then check for collisions. This can result in clusters overlapping.
    // Check clusters for collisions, and merge if needed.
    /*
    nclusters = clusters.size();
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
                --k;
            }
        }
    }
    */
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
    while (!completed) {
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
        }
        else {
            evolve(wld);
            ++step;
        }
        if (display) {
            i64 now = fenster_time();
            if (now - prev_render > 1000/FPS || completed) {
                prev_render = now;
                render(wld, &f);
            }
            //fenster_sleep(500);
        }
    }

    if (display) {
        if (completed) {
            // Display the final configuration.
            //fenster_loop(&f);
            // Pause before closing window.
            //fenster_sleep(1000);
        }
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
        return 1;
    }
    if (density < 0.0 || density >= 1.0) {
        fprintf(stderr, "Invalid density.\n");
        return 1;
    }
    if (threshold < 0.0 || threshold > 1.0) {
        fprintf(stderr, "Invalid threshold.\n");
        return 1;
    }


    srand(time(NULL));

    world wld;
    i32 min_steps = 0;
    i32 max_steps = 0;
    f64 mean_steps = 0.0;
    f64 sum_sqr_steps = 0.0;
    f64 stddev_steps = 0.0;
    for (i32 i = 1; i <= trials; ++i) {
        printf("Running trial %d/%d...\n", i, trials);
        populate(&wld, w, h, density);
        i32 steps = run(&wld, threshold, display, false);
        destroy(&wld);
        min_steps = (i == 1) ? steps : std::min(min_steps, steps);
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
