
#include "../rays/bvh.h"
#include "debug.h"
#include <stack>

#define N_BINS (16)

namespace PT {

// construct BVH hierarchy given a vector of prims
template<typename Primitive>
void BVH<Primitive>::build(std::vector<Primitive>&& prims, size_t max_leaf_size) {

    // NOTE (PathTracer):
    // This BVH is parameterized on the type of the primitive it contains. This allows
    // us to build a BVH over any type that defines a certain interface. Specifically,
    // we use this to both build a BVH over triangles within each Tri_Mesh, and over
    // a variety of Objects (which might be Tri_Meshes, Spheres, etc.) in Pathtracer.
    //
    // The Primitive interface must implement these two functions:
    //      BBox bbox() const;
    //      Trace hit(const Ray& ray) const;
    // Hence, you may call bbox() and hit() on any value of type Primitive.

    // Keep these two lines of code in your solution. They clear the list of nodes and
    // initialize member variable 'primitives' as a vector of the scene prims
    nodes.clear();
    primitives = std::move(prims);

    // TODO (PathTracer): Task 3
    // Modify the code ahead to construct a BVH from the given vector of primitives and maximum leaf
    // size configuration.
    //
    // Please use the SAH as described in class.  We recomment the binned build from lecture.
    // In general, here is a rough sketch:
    //
    //  For each axis X,Y,Z:
    //     Try possible splits along axis, evaluate SAH for each
    //  Take minimum cost across all axes.
    //  Partition primitives into a left and right child group
    //  Compute left and right child bboxes
    //  Make the left and right child nodes.
    //
    //
    // While a BVH is conceptually a tree structure, the BVH class uses a single vector (nodes)
    // to store all the nodes. Therefore, BVH nodes don't contain pointers to child nodes,
    // but rather the indices of the
    // child nodes in this array. Hence, to get the child of a node, you have to
    // look up the child index in this vector (e.g. nodes[node.l]). Similarly,
    // to create a new node, don't allocate one yourself - use BVH::new_node, which
    // returns the index of a newly added node.
    //
    // As an example of how to make nodes, the starter code below builds a BVH with a
    // root node that encloses all the primitives and its two descendants at Level 2.
    // For now, the split is hardcoded such that the first primitive is put in the left
    // child of the root, and all the other primitives are in the right child.
    // There are no further descendants.

    // edge case
    if(primitives.empty()) {
        return;
    }

    // compute bounding box for all primitives
    BBox bb;
    for(size_t i = 0; i < primitives.size(); ++i) {
        bb.enclose(primitives[i].bbox());
    }

    // set up root node (root BVH). Notice that it contains all primitives.
    size_t root_node_addr = new_node();
    root_idx = root_node_addr;
    Node& node = nodes[root_node_addr];
    node.bbox = bb;
    node.start = 0;
    node.size = primitives.size();

    build_subtree(root_node_addr, max_leaf_size);
}

/*
    We expect that both `n.l` and `n.r` are initialized to `0` by `new_node()`
*/
template<typename Primitive>
void BVH<Primitive>::build_subtree(size_t node_addr, size_t max_leaf_size) {
    Node n = nodes[node_addr];
    if (n.size <= max_leaf_size) {
        n.l = n.r = 0;
        return;
    }
    int best_axis = -1, best_split = -1; float min_cost = FLT_MAX; // [0, best_split] + [best_split + 1, N_BINS - 1]
    std::vector<int> best_bin_cnt; // bin `id` of `i`-the primitive
    BBox best_left, best_right;
    for (int axis = 0; axis < 3; ++axis) {
        float min = n.bbox.min[axis], max = n.bbox.max[axis];
        std::vector<int> cnt_bin; // number of elements in `i`-th bin
        std::vector<int> bin_cnt; // bin `id` of `i`-the primitive
        std::vector<BBox> bbox_bin;
        for (size_t i = 0; i < n.size; ++i)
            bin_cnt.push_back(0);
        for (size_t i = 0; i < N_BINS; ++i)
            cnt_bin.push_back(0), bbox_bin.push_back(BBox());
        for (size_t addr = n.start; addr < n.start + n.size; ++addr) {
            BBox b = primitives[addr].bbox();
            // Use center of bbox as centroid
            int bin = ((b.max[axis] + b.min[axis]) / 2 - min) * N_BINS / (max - min);
            bin = bin < 0 ? 0 : (bin >= N_BINS ? N_BINS - 1 : bin);
            bin_cnt[addr - n.start] = bin;
            cnt_bin[bin] += 1;
            bbox_bin[bin].enclose(b);
        }
        std::vector<BBox> left_bbox; // bbox of the first `i` bboxes
        std::vector<BBox> right_bbox; // bbox of the last `i` bboxes
        std::vector<int> left_sum; // total number of elements in the first `i` bboxes
        std::vector<int> right_sum; // total number of elements in the last `i` bboxes
        left_bbox.push_back(BBox()), right_bbox.push_back(BBox());
        left_sum.push_back(0), right_sum.push_back(0);
        for (size_t i = 0; i < N_BINS; ++i) {
            left_bbox.push_back(left_bbox[i]); left_bbox[i + 1].enclose(bbox_bin[i]);
            right_bbox.push_back(right_bbox[i]); right_bbox[i + 1].enclose(bbox_bin[N_BINS - i - 1]);
            left_sum.push_back(left_sum[i] + cnt_bin[i]);
            right_sum.push_back(right_sum[i] + cnt_bin[N_BINS - i - 1]);
        }
        for (size_t i = 1; i < N_BINS; ++i) {
            float cost = left_bbox[i].surface_area() * left_sum[i] + right_bbox[N_BINS - i].surface_area() * right_sum[N_BINS - i];
            if (cost < min_cost) {
                min_cost = cost;
                if (best_axis != axis)
                    best_bin_cnt = bin_cnt;
                best_axis = axis, best_split = i;
                best_left = left_bbox[i], best_right = right_bbox[N_BINS - i];
            }
        }
    }
    std::vector<Primitive> left, right;
    for (size_t i = 0; i < n.size; ++i) {
        if (best_bin_cnt[i] <= best_split)
            left.push_back(std::move(primitives[n.start + i]));
        else
            right.push_back(std::move(primitives[n.start + i]));
    }
    for (size_t i = 0; i < left.size(); ++i)
        primitives[n.start + i] = std::move(left[i]);
    for (size_t i = 0; i < right.size(); ++i)
        primitives[n.start + left.size() + i] = std::move(right[i]);

    if (left.size() == 0 || right.size() == 0) {
        n.l = n.r = 0;
        return;
    }

    // create child nodes
    size_t node_addr_l = new_node();
    size_t node_addr_r = new_node();
    nodes[node_addr].l = node_addr_l, nodes[node_addr].r = node_addr_r;

    nodes[node_addr_l].bbox = best_left;
    nodes[node_addr_l].start = n.start;
    nodes[node_addr_l].size = left.size();
    build_subtree(node_addr_l, max_leaf_size);

    nodes[node_addr_r].bbox = best_right;
    nodes[node_addr_r].start = n.start + left.size();
    nodes[node_addr_r].size = right.size();
    build_subtree(node_addr_r, max_leaf_size);
}

template<typename Primitive>
void BVH<Primitive>::hit_subtree(const Ray& ray, size_t node_addr, Trace& closest) const {
    Node n = nodes[node_addr];
    if (!(n.l && n.r)) {
        for (size_t idx = n.start; idx < n.start + n.size; ++idx) {
            const Primitive& prim = primitives[idx];
            Trace hit = prim.hit(ray);
            closest = Trace::min(closest, hit);
        }
    } else {
        Vec2 times_l = ray.dist_bounds, times_r = ray.dist_bounds;
        bool l_intersect = nodes[n.l].bbox.hit(ray, times_l);
        l_intersect = l_intersect && (!closest.hit || (times_l.x < closest.distance));
        bool r_intersect = nodes[n.r].bbox.hit(ray, times_r);
        r_intersect = r_intersect && (!closest.hit || (times_r.x < closest.distance));
        if (!l_intersect && !r_intersect)
            return;
        if (!l_intersect)
            hit_subtree(ray, n.r, closest);
        if (!r_intersect)
            hit_subtree(ray, n.l, closest);
        if (times_l.x < times_r.x) {
            hit_subtree(ray, n.l, closest);
            r_intersect = r_intersect && (!closest.hit || (times_r.x < closest.distance));
            if (r_intersect)
                hit_subtree(ray, n.r, closest);
        }
        else {
            hit_subtree(ray, n.r, closest);
            l_intersect = l_intersect && (!closest.hit || (times_l.x < closest.distance));
            if (l_intersect)
                hit_subtree(ray, n.l, closest);
        }
    }
}

template<typename Primitive>
Trace BVH<Primitive>::hit(const Ray& ray) const {

    // TODO (PathTracer): Task 3
    // Implement ray - BVH intersection test. A ray intersects
    // with a BVH aggregate if and only if it intersects a primitive in
    // the BVH that is not an aggregate.

    // The starter code simply iterates through all the primitives.
    // Again, remember you can use hit() on any Primitive value.

    Trace closest;
    if (nodes.empty())
        return closest;
    Vec2 times = ray.dist_bounds;
    if (nodes[root_idx].bbox.hit(ray, times))
        hit_subtree(ray, root_idx, closest);
    return closest;
}

template<typename Primitive>
BVH<Primitive>::BVH(std::vector<Primitive>&& prims, size_t max_leaf_size) {
    build(std::move(prims), max_leaf_size);
}

template<typename Primitive>
BVH<Primitive> BVH<Primitive>::copy() const {
    BVH<Primitive> ret;
    ret.nodes = nodes;
    ret.primitives = primitives;
    ret.root_idx = root_idx;
    return ret;
}

template<typename Primitive>
bool BVH<Primitive>::Node::is_leaf() const {
    return l == r;
}

template<typename Primitive>
size_t BVH<Primitive>::new_node(BBox box, size_t start, size_t size, size_t l, size_t r) {
    Node n;
    n.bbox = box;
    n.start = start;
    n.size = size;
    n.l = l;
    n.r = r;
    nodes.push_back(n);
    return nodes.size() - 1;
}

template<typename Primitive>
BBox BVH<Primitive>::bbox() const {
    return nodes[root_idx].bbox;
}

template<typename Primitive>
std::vector<Primitive> BVH<Primitive>::destructure() {
    nodes.clear();
    return std::move(primitives);
}

template<typename Primitive>
void BVH<Primitive>::clear() {
    nodes.clear();
    primitives.clear();
}

template<typename Primitive>
size_t BVH<Primitive>::visualize(GL::Lines& lines, GL::Lines& active, size_t level,
                                 const Mat4& trans) const {

    std::stack<std::pair<size_t, size_t>> tstack;
    tstack.push({root_idx, 0});
    size_t max_level = 0;

    if(nodes.empty()) return max_level;

    while(!tstack.empty()) {

        auto [idx, lvl] = tstack.top();
        max_level = std::max(max_level, lvl);
        const Node& node = nodes[idx];
        tstack.pop();

        Vec3 color = lvl == level ? Vec3(1.0f, 0.0f, 0.0f) : Vec3(1.0f);
        GL::Lines& add = lvl == level ? active : lines;

        BBox box = node.bbox;
        box.transform(trans);
        Vec3 min = box.min, max = box.max;

        auto edge = [&](Vec3 a, Vec3 b) { add.add(a, b, color); };

        edge(min, Vec3{max.x, min.y, min.z});
        edge(min, Vec3{min.x, max.y, min.z});
        edge(min, Vec3{min.x, min.y, max.z});
        edge(max, Vec3{min.x, max.y, max.z});
        edge(max, Vec3{max.x, min.y, max.z});
        edge(max, Vec3{max.x, max.y, min.z});
        edge(Vec3{min.x, max.y, min.z}, Vec3{max.x, max.y, min.z});
        edge(Vec3{min.x, max.y, min.z}, Vec3{min.x, max.y, max.z});
        edge(Vec3{min.x, min.y, max.z}, Vec3{max.x, min.y, max.z});
        edge(Vec3{min.x, min.y, max.z}, Vec3{min.x, max.y, max.z});
        edge(Vec3{max.x, min.y, min.z}, Vec3{max.x, max.y, min.z});
        edge(Vec3{max.x, min.y, min.z}, Vec3{max.x, min.y, max.z});

        if(node.l && node.r) {
            tstack.push({node.l, lvl + 1});
            tstack.push({node.r, lvl + 1});
        } else {
            for(size_t i = node.start; i < node.start + node.size; i++) {
                size_t c = primitives[i].visualize(lines, active, level - lvl, trans);
                max_level = std::max(c, max_level);
            }
        }
    }
    return max_level;
}

} // namespace PT
