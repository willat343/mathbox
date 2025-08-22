#ifndef MATHBOX_IMPL_TRANSFORM_MANAGER_HPP
#define MATHBOX_IMPL_TRANSFORM_MANAGER_HPP

#include <cppbox/exceptions.hpp>

#include "mathbox/transform_manager.hpp"

namespace math {

template<int D_>
    requires(math::is_2d_or_3d<D_>)
Transform<D_>::Transform(const std::string& parent_frame_, const std::string& child_frame_, const Pose<D>& transform_)
    : parent_frame_(parent_frame_),
      child_frame_(child_frame_),
      transform_(transform_) {
    throw_if(parent_frame_ == child_frame_ && !transform_.matrix().isIdentity(),
            "Transform has parent_frame == child_frame but is not identity.");
}

template<int D_>
    requires(math::is_2d_or_3d<D_>)
const std::string& Transform<D_>::parent_frame() const {
    return parent_frame_;
}

template<int D_>
    requires(math::is_2d_or_3d<D_>)
const std::string& Transform<D_>::child_frame() const {
    return child_frame_;
}

template<int D_>
    requires(math::is_2d_or_3d<D_>)
auto Transform<D_>::transform() const -> const Pose<D>& {
    return transform_;
}

template<int D_>
    requires(math::is_2d_or_3d<D_>)
void TransformManager<D_>::add_transform(const std::string& parent_frame, const std::string& child_frame,
        const Pose<D>& transform) {
    // Throw error if the parent_frame and child_frame are the same:
    throw_if(parent_frame == child_frame, "TransformManager cannot add transform where parent frame (" + parent_frame +
                                                  ") and child frame (" + child_frame + ") are the same.");
    // Throw error if child_frame exists and is not a root node
    throw_if(has_child_frame(child_frame),
            "TransformManager cannot add " + parent_frame + "->" + child_frame + " transform because " + child_frame +
                    " already exists with parent " + get_frame(child_frame).parent_transform()->parent_frame().name() +
                    ". Child frames may only have one parent.");
    // Get or create parent_frame with book-keeping for the new child_frame
    FrameNode& frame = get_or_create_frame_for_child_frame(parent_frame, child_frame);
    // Add the transform to the parent FrameNode
    frame.add_transform(child_frame, transform);
    // Move existing tree (with root frame == child_frame) under the new child FrameNode
    for (auto it = roots().begin(); it != roots().end(); ++it) {
        if (it->name() == child_frame) {
            // Move root FrameNode to child FrameNode (pointers and references remain valid)
            frame.transform(child_frame).child_frame() = std::move(*it);
            // Set parent TransformEdge of child FrameNode
            frame.transform(child_frame).child_frame().parent_transform() = &frame.transform(child_frame);
            // Update book-keeping for each FrameNode above frame
            TransformEdge* current_transform = frame.parent_transform();
            while (current_transform != nullptr) {
                current_transform->parent_frame().add_next_child_frame_for_frame(child_frame,
                        current_transform->child_frame().name());
                current_transform = current_transform->parent_frame().parent_transform();
            }
            // Delete moved root FrameNode (invalidating iterators)
            roots().erase(it);
            // Terminate loop as can only be one root FrameNode to move (and iterator is invalidated by erase call)
            break;
        }
    }
}

template<int D_>
    requires(math::is_2d_or_3d<D_>)
inline void TransformManager<D_>::add_transform(const Transform<D>& transform) {
    return add_transform(transform.parent_frame(), transform.child_frame(), transform.transform());
}

template<int D_>
    requires(math::is_2d_or_3d<D_>)
bool TransformManager<D_>::has_child_frame(const std::string& frame) const {
    for (auto it = roots().cbegin(); it != roots().cend(); ++it) {
        if (it->has_child_frame(frame)) {
            return true;
        }
    }
    return false;
}

template<int D_>
    requires(math::is_2d_or_3d<D_>)
bool TransformManager<D_>::has_frame(const std::string& frame) const {
    for (auto it = roots().cbegin(); it != roots().cend(); ++it) {
        if (it->has_frame(frame)) {
            return true;
        }
    }
    return false;
}

template<int D_>
    requires(math::is_2d_or_3d<D_>)
bool TransformManager<D_>::has_transform(const std::string& parent_frame, const std::string& child_frame) const {
    for (auto it = roots().cbegin(); it != roots().cend(); ++it) {
        if (it->has_frame(parent_frame) && it->has_frame(child_frame)) {
            return true;
        }
    }
    return false;
}

template<int D_>
    requires(math::is_2d_or_3d<D_>)
inline std::size_t TransformManager<D_>::num_roots() const {
    return roots().size();
}

template<int D_>
    requires(math::is_2d_or_3d<D_>)
std::deque<std::string> TransformManager<D_>::root_frames() const {
    const std::size_t size = num_roots();
    std::deque<std::string> root_names_(size);
    for (std::size_t i = 0; i < size; ++i) {
        root_names_[i] = roots()[i].name();
    }
    return root_names_;
}

template<int D_>
    requires(math::is_2d_or_3d<D_>)
auto TransformManager<D_>::transform(const std::string& parent_frame, const std::string& child_frame) const -> Pose<D> {
    // Get root FrameNode, which checks for existence of parent FrameNode
    const FrameNode* current_frame = &get_root_for_frame(parent_frame);
    // Catch early exit if frames are identical
    if (parent_frame == child_frame) {
        return Pose<D>::Identity();
    }
    // Check for existence of child_frame
    throw_if(!current_frame->has_frame(child_frame), "TransformManager does not have frame " + child_frame +
                                                             " in tree rooted at frame " + current_frame->name() +
                                                             " where frame " + parent_frame + " was found.");
    // Compute transforms from last common parent FrameNode to the parent FrameNode and to the child FrameNode
    bool found_common{false};
    Pose<D> transform_common_parent = Pose<D>::Identity();
    Pose<D> transform_common_child = Pose<D>::Identity();
    while (!found_common) {
        if (current_frame->name() == parent_frame) {
            // Parent FrameNode found first
            found_common = true;
            // Acquire transform from common/parent FrameNode to child frame
            transform_common_child = current_frame->next_transform(child_frame).transform_to(child_frame);
        } else if (current_frame->name() == child_frame) {
            // Child FrameNode found first
            found_common = true;
            // Acquire transform from common/child FrameNode to parent frame
            transform_common_parent = current_frame->next_transform(parent_frame).transform_to(parent_frame);
        } else {
            // Niether parent nor child found so continue to traverse to find last common parent FrameNode
            const std::string next_child_for_parent_frame = current_frame->next_child_frame_name(parent_frame);
            const std::string next_child_for_child_frame = current_frame->next_child_frame_name(child_frame);
            if (next_child_for_parent_frame == next_child_for_child_frame) {
                // Update current frame to continue searching for last common frame
                current_frame = &current_frame->transform(next_child_for_parent_frame).child_frame();
            } else {
                // Common parent FrameNode found but is neither the parent nor child frame
                found_common = true;
                // Acquire transform from common FrameNode to parent frame
                transform_common_parent =
                        current_frame->transform(next_child_for_parent_frame).transform_to(parent_frame);
                // Acquire transform from common FrameNode to child frame
                transform_common_child = current_frame->transform(next_child_for_child_frame).transform_to(child_frame);
            }
        }
    }
    // Compute transform
    return transform_common_parent.inverse() * transform_common_child;
}

template<int D_>
    requires(math::is_2d_or_3d<D_>)
inline auto TransformManager<D_>::transforms() const -> std::deque<Transform<D>> {
    std::deque<Transform<D>> transforms_;
    for (const FrameNode& root : roots()) {
        root.append_transforms(transforms_);
    }
    return transforms_;
}

template<int D_>
    requires(math::is_2d_or_3d<D_>)
TransformManager<D_>::FrameNode::FrameNode(const std::string& frame_name_, TransformEdge* parent_transform_)
    : frame_name_(frame_name_),
      parent_transform_(parent_transform_) {}

template<int D_>
    requires(math::is_2d_or_3d<D_>)
inline void TransformManager<D_>::FrameNode::add_next_child_frame_for_frame(const std::string& frame_name,
        const std::string& child_frame_name) {
    next_child_frame_for_frame.emplace(frame_name, child_frame_name);
}

template<int D_>
    requires(math::is_2d_or_3d<D_>)
inline void TransformManager<D_>::FrameNode::add_transform(const std::string& child_frame_name,
        const Pose<D>& transform) {
    add_next_child_frame_for_frame(child_frame_name, child_frame_name);
    // By using dynamic memory, the pointer remains valid after emplacing into the map
    child_transforms_.emplace(child_frame_name, std::make_unique<TransformEdge>(*this, child_frame_name, transform));
}

template<int D_>
    requires(math::is_2d_or_3d<D_>)
void TransformManager<D_>::FrameNode::append_transforms(std::deque<Transform<D>>& transforms_) const {
    for (const auto& [child_frame_name, child_transform] : child_transforms()) {
        transforms_.emplace_back(name(), child_frame_name, child_transform->transform());
        child_transform->child_frame().append_transforms(transforms_);
    }
}

template<int D_>
    requires(math::is_2d_or_3d<D_>)
inline std::deque<std::string> TransformManager<D_>::FrameNode::child_frame_names() const {
    std::deque<std::string> child_frame_names_;
    for (auto it = child_transforms_.cbegin(); it != child_transforms_.cend(); ++it) {
        child_frame_names_.emplace_back(it->first);
    }
    return child_frame_names_;
}

template<int D_>
    requires(math::is_2d_or_3d<D_>)
inline auto TransformManager<D_>::FrameNode::child_transforms() const
        -> const std::map<std::string, std::unique_ptr<TransformEdge>>& {
    return child_transforms_;
}

template<int D_>
    requires(math::is_2d_or_3d<D_>)
inline auto TransformManager<D_>::FrameNode::child_transforms()
        -> std::map<std::string, std::unique_ptr<TransformEdge>>& {
    return const_cast<std::map<std::string, std::unique_ptr<TransformEdge>>&>(std::as_const(*this).child_transforms());
}

template<int D_>
    requires(math::is_2d_or_3d<D_>)
inline const std::string& TransformManager<D_>::FrameNode::name() const {
    return frame_name_;
}

template<int D_>
    requires(math::is_2d_or_3d<D_>)
inline bool TransformManager<D_>::FrameNode::has_frame(const std::string& query_frame_name) const {
    return query_frame_name == name() || has_child_frame(query_frame_name);
}

template<int D_>
    requires(math::is_2d_or_3d<D_>)
inline bool TransformManager<D_>::FrameNode::has_child_frame(const std::string& query_frame_name) const {
    return next_child_frame_for_frame.find(query_frame_name) != next_child_frame_for_frame.end();
}

template<int D_>
    requires(math::is_2d_or_3d<D_>)
inline const std::string& TransformManager<D_>::FrameNode::next_child_frame_name(
        const std::string& query_frame_name) const {
    return next_child_frame_for_frame.at(query_frame_name);
}

template<int D_>
    requires(math::is_2d_or_3d<D_>)
inline const typename TransformManager<D_>::TransformEdge& TransformManager<D_>::FrameNode::next_transform(
        const std::string& query_frame_name) const {
    return *child_transforms_.at(next_child_frame_name(query_frame_name));
}

template<int D_>
    requires(math::is_2d_or_3d<D_>)
inline const typename TransformManager<D_>::TransformEdge* const& TransformManager<D_>::FrameNode::parent_transform()
        const {
    return parent_transform_;
}

template<int D_>
    requires(math::is_2d_or_3d<D_>)
inline typename TransformManager<D_>::TransformEdge*& TransformManager<D_>::FrameNode::parent_transform() {
    return const_cast<TransformEdge*&>(std::as_const(*this).parent_transform());
}

template<int D_>
    requires(math::is_2d_or_3d<D_>)
inline const typename TransformManager<D_>::TransformEdge& TransformManager<D_>::FrameNode::transform(
        const std::string& child_frame_name) const {
    return *child_transforms_.at(child_frame_name);
}

template<int D_>
    requires(math::is_2d_or_3d<D_>)
inline typename TransformManager<D_>::TransformEdge& TransformManager<D_>::FrameNode::transform(
        const std::string& child_frame_name) {
    return const_cast<TransformEdge&>(std::as_const(*this).transform(child_frame_name));
}

template<int D_>
    requires(math::is_2d_or_3d<D_>)
TransformManager<D_>::TransformEdge::TransformEdge(const FrameNode& parent_frame_, const std::string& child_frame_name,
        const Pose<D>& transform_)
    : parent_frame_(parent_frame_),
      child_frame_(FrameNode{child_frame_name, this}),
      transform_(transform_) {}

template<int D_>
    requires(math::is_2d_or_3d<D_>)
inline const typename TransformManager<D_>::FrameNode& TransformManager<D_>::TransformEdge::child_frame() const {
    return child_frame_;
}

template<int D_>
    requires(math::is_2d_or_3d<D_>)
inline typename TransformManager<D_>::FrameNode& TransformManager<D_>::TransformEdge::child_frame() {
    return const_cast<FrameNode&>(std::as_const(*this).child_frame());
}

template<int D_>
    requires(math::is_2d_or_3d<D_>)
inline const typename TransformManager<D_>::FrameNode& TransformManager<D_>::TransformEdge::parent_frame() const {
    return parent_frame_;
}

template<int D_>
    requires(math::is_2d_or_3d<D_>)
inline typename TransformManager<D_>::FrameNode& TransformManager<D_>::TransformEdge::parent_frame() {
    return const_cast<FrameNode&>(std::as_const(*this).parent_frame());
}

template<int D_>
    requires(math::is_2d_or_3d<D_>)
inline auto TransformManager<D_>::TransformEdge::transform() const -> const Pose<D>& {
    return transform_;
}

template<int D_>
    requires(math::is_2d_or_3d<D_>)
auto TransformManager<D_>::TransformEdge::transform_to(const std::string& child_frame_name) const -> Pose<D> {
    const TransformEdge* current_transform = this;
    Pose<D> transform_to_child = current_transform->transform();
    while (current_transform->child_frame().name() != child_frame_name) {
        current_transform = &current_transform->child_frame().next_transform(child_frame_name);
        transform_to_child = transform_to_child * current_transform->transform();
    }
    return transform_to_child;
}

template<int D_>
    requires(math::is_2d_or_3d<D_>)
inline typename TransformManager<D_>::FrameNode& TransformManager<D_>::create_root(const std::string& frame_name) {
    return roots_.emplace_back(frame_name, nullptr);
}

template<int D_>
    requires(math::is_2d_or_3d<D_>)
inline const TransformManager<D_>::FrameNode& TransformManager<D_>::get_frame(const std::string& frame_name) const {
    // Check if frame exists
    throw_if(!has_frame(frame_name), "TransformManager does not have Frame " + frame_name + ".");
    // Find existing parent FrameNode
    const FrameNode* current_frame = &get_root_for_frame(frame_name);
    while (current_frame->name() != frame_name) {
        current_frame = &current_frame->transform(current_frame->next_child_frame_name(frame_name)).child_frame();
    }
    return *current_frame;
}

template<int D_>
    requires(math::is_2d_or_3d<D_>)
const typename TransformManager<D_>::FrameNode& TransformManager<D_>::get_root_for_frame(
        const std::string& frame_name) const {
    const std::size_t num_roots_ = num_roots();
    for (std::size_t i = 0; i < num_roots_; ++i) {
        if (roots()[i].has_frame(frame_name)) {
            return roots()[i];
        }
    }
    throw_here("TransformManager does not have Frame " + frame_name + ".");
}

template<int D_>
    requires(math::is_2d_or_3d<D_>)
inline typename TransformManager<D_>::FrameNode& TransformManager<D_>::get_root_for_frame(
        const std::string& frame_name) {
    return const_cast<FrameNode&>(std::as_const(*this).get_root_for_frame(frame_name));
}

template<int D_>
    requires(math::is_2d_or_3d<D_>)
typename TransformManager<D_>::FrameNode& TransformManager<D_>::get_or_create_frame_for_child_frame(
        const std::string& frame_name, const std::string& new_child_frame_name) {
    // Check if frame exists
    if (has_frame(frame_name)) {
        // Find existing parent FrameNode
        FrameNode* current_frame = &get_root_for_frame(frame_name);
        while (current_frame->name() != frame_name) {
            const std::string& next_child_frame = current_frame->next_child_frame_name(frame_name);
            // Add new_child_frame_name to each FrameNode's next_child_frame_for_frame map
            current_frame->add_next_child_frame_for_frame(new_child_frame_name, next_child_frame);
            current_frame = &current_frame->transform(next_child_frame).child_frame();
        }
        return *current_frame;
    } else {
        // Create new tree with transform
        return create_root(frame_name);
    }
}

template<int D_>
    requires(math::is_2d_or_3d<D_>)
inline const std::deque<typename TransformManager<D_>::FrameNode>& TransformManager<D_>::roots() const {
    return roots_;
}

template<int D_>
    requires(math::is_2d_or_3d<D_>)
inline std::deque<typename TransformManager<D_>::FrameNode>& TransformManager<D_>::roots() {
    return const_cast<std::deque<typename TransformManager<D_>::FrameNode>&>(std::as_const(*this).roots());
}

}

#if !MATHBOX_HEADER_ONLY
namespace math {

extern template class TransformManager<2>;
extern template class TransformManager<3>;

}
#endif

#endif
