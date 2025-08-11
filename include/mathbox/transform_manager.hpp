#ifndef MATHBOX_TRANSFORM_MANAGER_HPP
#define MATHBOX_TRANSFORM_MANAGER_HPP

#include <Eigen/Geometry>
#include <deque>
#include <map>
#include <memory>
#include <string>
#include <type_traits>

#include "mathbox/traits.hpp"
#include "mathbox/types.hpp"

namespace math {

/**
 * @brief Transform class.
 *
 * @tparam D_ dimension
 */
template<int D_>
    requires(math::is_2d_or_3d<D_>)
class Transform {
public:
    /**
     * @brief Dimension D (either 2 for 2D or 3 for 3D).
     *
     */
    static constexpr int D = D_;

    explicit Transform(const std::string& parent_frame_, const std::string& child_frame_, const Pose<D>& transform_);

    const std::string& parent_frame() const;

    const std::string& child_frame() const;

    const Pose<D>& transform() const;

private:
    std::string parent_frame_;
    std::string child_frame_;
    Pose<D> transform_;
};

template<int D_>
    requires(math::is_2d_or_3d<D_>)
class TransformManager {
public:
    /**
     * @brief Dimension D (either 2 for 2D or 3 for 3D).
     *
     */
    static constexpr int D = D_;

    /**
     * @brief Add a transform from parent_frame to child_frame, T_P_C, i.e. the pose of the child frame in the parent
     * frame.
     *
     * @param parent_frame
     * @param child_frame
     * @param transform
     */
    void add_transform(const std::string& parent_frame, const std::string& child_frame, const Pose<D>& transform);

    /**
     * @brief Convenience overload.
     *
     * @param transform
     */
    void add_transform(const Transform<D>& transform);

    /**
     * @brief Check if frame exists as a child frame in a transform tree (i.e. not a root FrameNode).
     *
     * @param frame
     * @return true
     * @return false
     */
    bool has_child_frame(const std::string& frame) const;

    /**
     * @brief Check if frame exists.
     *
     * @param frame
     * @return true
     * @return false
     */
    bool has_frame(const std::string& frame) const;

    /**
     * @brief Check if transform between two frames exists. They only need to be part of the same connected tree for
     * this to be true (they do not need to be directly connected, and the parent-child order does not matter).
     *
     * @param parent_frame
     * @param child_frame
     * @return true
     * @return false
     */
    bool has_transform(const std::string& parent_frame, const std::string& child_frame) const;

    /**
     * @brief Get the number of root Frames (i.e. number of disconnected transform trees).
     *
     * @return std::size_t
     */
    std::size_t num_roots() const;

    /**
     * @brief Get the root frame names.
     *
     * @return std::deque<std::string>
     */
    std::deque<std::string> root_frames() const;

    /**
     * @brief Get the transform from parent_frame to child_frame, T_P_C, i.e. the pose of the child frame in the parent
     * frame.
     *
     * @param parent_frame
     * @param child_frame
     * @return Pose<D>
     */
    Pose<D> transform(const std::string& parent_frame, const std::string& child_frame) const;

    /**
     * @brief Get all transform pairs in the trees.
     *
     * @return std::deque<Transform<D>>
     */
    std::deque<Transform<D>> transforms() const;

private:
    // Forward declaration
    class TransformEdge;

    /**
     * @brief A FrameNode is defined by a name, and can be considered a node in the transform tree. It may contain a
     * number of Transforms to child FrameNodes. All FrameNode names below this FrameNode in the transform tree are
     * known to the FrameNode.
     *
     */
    class FrameNode {
    public:
        /**
         * @brief Construct a new FrameNode.
         *
         * @param frame_name_
         */
        explicit FrameNode(const std::string& frame_name_, TransformEdge* parent_transform_);

        /**
         * @brief Add a frame->next_child_frame mapping to FrameNode so that frame can be found during a query.
         *
         * @param frame_name
         * @param child_frame_name
         */
        void add_next_child_frame_for_frame(const std::string& frame_name, const std::string& child_frame_name);

        /**
         * @brief Add a transform from this frame to the child_frame.
         *
         * @param child_frame_name
         * @param transform
         */
        void add_transform(const std::string& child_frame_name, const Pose<D>& transform);

        /**
         * @brief Append all child transforms recursively to the container.
         *
         * @param transforms_
         */
        void append_transforms(std::deque<Transform<D>>& transforms_) const;

        /**
         * @brief Get the child frame names.
         *
         * @return std::deque<std::string>
         */
        std::deque<std::string> child_frame_names() const;

        /**
         * @brief Get the child transforms.
         *
         * @return std::map<std::string, std::unique_ptr<TransformEdge>>&
         */
        const std::map<std::string, std::unique_ptr<TransformEdge>>& child_transforms() const;

        /**
         * @brief Get the mutable child transforms.
         *
         * @return std::map<std::string, std::unique_ptr<TransformEdge>>&
         */
        std::map<std::string, std::unique_ptr<TransformEdge>>& child_transforms();

        /**
         * @brief Get frame name.
         *
         * @return const std::string&
         */
        const std::string& name() const;

        /**
         * @brief Check if query_frame_name is a child frame.
         *
         * @param query_frame_name
         * @return true
         * @return false
         */
        bool has_child_frame(const std::string& query_frame_name) const;

        /**
         * @brief Check if query_frame_name is this frame or a child frame.
         *
         * @param query_frame_name
         * @return true
         * @return false
         */
        bool has_frame(const std::string& query_frame_name) const;

        /**
         * @brief Get the next child frame name that would be needed to get to query_frame_name during tree traversal.
         *
         * @param query_frame_name
         * @return const std::string&
         */
        const std::string& next_child_frame_name(const std::string& query_frame_name) const;

        /**
         * @brief Directly get the next TransformEdge for the traversal to query_frame_name.
         *
         * @param query_frame_name
         * @return const TransformEdge&
         */
        const TransformEdge& next_transform(const std::string& query_frame_name) const;

        /**
         * @brief Get the parent TransformEdge.
         *
         * @return const TransformEdge* const&
         */
        const TransformEdge* const& parent_transform() const;

        /**
         * @brief Get the parent TransformEdge.
         *
         * @return TransformEdge*&
         */
        TransformEdge*& parent_transform();

        /**
         * @brief Get the TransformEdge for child_frame_name.
         *
         * @param child_frame_name
         * @return const TransformEdge&
         */
        const TransformEdge& transform(const std::string& child_frame_name) const;

        /**
         * @brief Get the TransformEdge for child_frame_name.
         *
         * @param child_frame_name
         * @return TransformEdge&
         */
        TransformEdge& transform(const std::string& child_frame_name);

    private:
        /**
         * @brief FrameNode name.
         *
         */
        std::string frame_name_;

        /**
         * @brief Parent transform to this frame.
         *
         */
        TransformEdge* parent_transform_;

        /**
         * @brief Direct transforms to child frames from this frame.
         *
         */
        std::map<std::string, std::unique_ptr<TransformEdge>> child_transforms_;

        /**
         * @brief Map of the name of the child frame being queried to the name of the child frame from this node where
         * it can be found. This map includes all the child frames.
         *
         * This avoids the need for tree search algorithms.
         *
         */
        std::map<std::string, std::string> next_child_frame_for_frame;
    };

    /**
     * @brief A TransformEdge is defined by a Pose, and may be considered an edge in the transform tree.
     *
     */
    class TransformEdge {
    public:
        /**
         * @brief Construct a new TransformEdge.
         *
         * @param parent_frame_
         * @param child_frame_
         * @param transform_
         */
        explicit TransformEdge(const FrameNode& parent_frame_, const std::string& child_frame_name,
                const Pose<D>& transform_);

        /**
         * @brief Get the child FrameNode.
         *
         * @return const FrameNode&
         */
        const FrameNode& child_frame() const;

        /**
         * @brief Get the child FrameNode.
         *
         * @return FrameNode&
         */
        FrameNode& child_frame();

        /**
         * @brief Get the parent FrameNode.
         *
         * @return const FrameNode&
         */
        const FrameNode& parent_frame() const;

        /**
         * @brief Get te parent FrameNode.
         *
         * @return FrameNode&
         */
        FrameNode& parent_frame();

        /**
         * @brief Get the transform.
         *
         * @return const Pose<D>&
         */
        const Pose<D>& transform() const;

        /**
         * @brief Compute the transform from this transform (i.e. including this transform from the parent FrameNode) to
         * a child frame that exists.
         *
         * @param child_frame_name
         * @return Pose<D>
         */
        Pose<D> transform_to(const std::string& child_frame_name) const;

    private:
        /**
         * @brief Parent frame.
         *
         */
        const FrameNode& parent_frame_;

        /**
         * @brief Child frame.
         *
         */
        FrameNode child_frame_;

        /**
         * @brief TransformEdge from parent frame to child frame, T_P_C, i.e. pose of child frame in parent frame.
         *
         */
        Pose<D> transform_;
    };

private:
    /**
     * @brief Create a root frame.
     *
     * @param frame
     */
    FrameNode& create_root(const std::string& frame_name);

    /**
     * @brief Get a frame.
     *
     * @param frame_name
     * @return const FrameNode&
     */
    const FrameNode& get_frame(const std::string& frame_name) const;

    /**
     * @brief Get the root FrameNode that contains frame.
     *
     * @param frame_name
     * @return const FrameNode&
     */
    const FrameNode& get_root_for_frame(const std::string& frame_name) const;

    /**
     * @brief Get the root FrameNode that contains frame.
     *
     * @param frame_name
     * @return FrameNode&
     */
    FrameNode& get_root_for_frame(const std::string& frame_name);

    /**
     * @brief Get or create frame object for a new transform to child_frame, recursively adding to next_child_for_frame
     * mappings (if necessary). The new transform is not created.
     *
     * @param frame
     * @return FrameNode&
     */
    FrameNode& get_or_create_frame_for_child_frame(const std::string& frame_name,
            const std::string& new_child_frame_name);

    /**
     * @brief Get the root Frames.
     *
     * @return const std::deque<FrameNode>&
     */
    const std::deque<FrameNode>& roots() const;

    /**
     * @brief Get the root Frames.
     *
     * @return std::deque<FrameNode>&
     */
    std::deque<FrameNode>& roots();

    /**
     * @brief Root nodes of the transform frame tree.
     *
     */
    std::deque<FrameNode> roots_;
};

}

#include "mathbox/impl/transform_manager.hpp"

#endif
