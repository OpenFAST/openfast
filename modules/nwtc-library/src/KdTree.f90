!> @module kdtree_module
!> @brief Provides a k-d tree implementation for efficient nearest neighbor searches.
!>
!> This module implements a static k-d tree with variable dimensionality.
!> The tree is built in-place from a given set of points, which is a highly 
!> efficient method in terms of memory and speed. It is ideal for scenarios 
!> where the set of points does not change after the initial construction.
module KdTree

use NWTC_Base
use NWTC_Library_Types

implicit none
private

public :: kdtree_build
public :: kdtree_destroy
public :: kdtree_nearest_neighbor
public :: kdtree_points_in_radius
public :: kdtree_update
public :: kdtree_sort_indices

contains

!> @brief Builds a k-d tree from a set of points.
!> @param[out] tree The k-d tree object to be constructed.
!> @param[in] points_in A (NDims, N) array of N points with NDims dimensions.
!> @param[in] n_active Optional: Number of active points to use (default: all points).
!> @param[in] n_max Optional: Maximum capacity for future updates (default: n_active).
subroutine kdtree_build(tree, points_in, n_max)
   type(KdTreeType), intent(OUT) :: tree
   real(ReKi), dimension(:, :), intent(IN) :: points_in
   integer(IntKi), intent(IN), optional :: n_max
   integer(IntKi) :: i, n_pts, n_capacity, n_dims

   ! Get dimensionality from input array
   n_dims = size(points_in, 1)
   n_pts = size(points_in, 2)
   
   ! Determine number of points in tree
   tree%NNodes = n_pts
   
   ! Determine maximum capacity
   if (present(n_max)) then
      n_capacity = max(n_max, tree%NNodes)
   else
      n_capacity = tree%NNodes
   end if

   allocate (tree%Nodes(n_dims, 0:n_capacity))
   tree%Nodes = 0.0_ReKi
   allocate (tree%OriginalIndices(n_capacity))
   tree%OriginalIndices(tree%NNodes:) = 0
   allocate (tree%Dimensions(n_capacity))
   tree%Dimensions(tree%NNodes:) = 0

   tree%Nodes(:, 1:tree%NNodes) = points_in(:, 1:tree%NNodes)
   tree%OriginalIndices(1:tree%NNodes) = [(i, i=1, tree%NNodes)]

   ! Start the recursive build process on the active nodes.
   call build_recursive(tree, 1, tree%NNodes, 0)
end subroutine kdtree_build

!> @brief Deallocates the memory used by a k-d tree object.
subroutine kdtree_destroy(tree)
   type(KdTreeType), intent(INOUT) :: tree
   if (ALLOCATED(tree%Nodes)) deallocate (tree%Nodes)
   if (ALLOCATED(tree%OriginalIndices)) deallocate (tree%OriginalIndices)
   if (ALLOCATED(tree%Dimensions)) deallocate (tree%Dimensions)
   tree%NNodes = 0
end subroutine kdtree_destroy

!> @brief Finds the index of the nearest neighbor to a query point.
!> @param[in] tree The constructed k-d tree.
!> @param[in] query_point The point to search for (must match tree dimensionality).
!> @return The original index of the nearest neighbor point.
function kdtree_nearest_neighbor(tree, query_point) result(nn_idx)
   type(KdTreeType), intent(IN) :: tree
   real(ReKi), dimension(:), intent(IN) :: query_point
   integer(IntKi) :: nn_idx
   real(ReKi) :: best_dist_sq

   if (tree%NNodes <= 0) then
      nn_idx = 0 ! Or some other indicator of error/empty tree
      return
   end if

   best_dist_sq = HUGE(1.0_ReKi)
   nn_idx = -1

   ! Start the recursive search. The root of the full tree is at the median index.
   call search_recursive(tree, 1, tree%NNodes, query_point, best_dist_sq, nn_idx)
end function kdtree_nearest_neighbor

!> @brief Finds all points within a given radius of a query point.
!> @param[in] tree The constructed k-d tree.
!> @param[in] query_point The point to search around (must match tree dimensionality).
!> @param[in] radius The search radius.
!> @param[out] indices Array of original indices of points within radius.
!> @param[out] n_found The number of points found within radius.
subroutine kdtree_points_in_radius(tree, query_point, radius, indices, n_found)
   type(KdTreeType), intent(IN) :: tree
   real(ReKi), dimension(:), intent(IN) :: query_point
   real(ReKi), intent(IN) :: radius
   integer, dimension(:), intent(OUT) :: indices
   integer, intent(OUT) :: n_found
   real(ReKi) :: radius_sq

   ! Indices must be capable of storing indices for all nodes
   if (size(indices) < tree%NNodes) then
      n_found = -1
      return
   end if

   ! If the tree is empty, return no nodes found
   if (tree%NNodes <= 0) then
      n_found = 0
      return
   end if

   n_found = 0
   radius_sq = radius**2

   ! Start the recursive search
   call radius_search_recursive(tree, 1, tree%NNodes, query_point, radius_sq, indices, n_found)
end subroutine kdtree_points_in_radius

!> @brief Sorts an array of indices in ascending order.
!> @param[inout] indices Array of indices to be sorted.
!> @details Uses insertion sort for small arrays (< 20 elements) and quicksort for larger arrays.
subroutine kdtree_sort_indices(indices)
   integer, dimension(:), intent(INOUT) :: indices
   integer :: n
   
   n = size(indices)
   
   if (n <= 1) return
   
   ! Use insertion sort for small arrays, quicksort for larger ones
   if (n < 20) then
      call insertion_sort(indices, 1, n)
   else
      call quicksort_indices(indices, 1, n)
   end if
end subroutine kdtree_sort_indices

!> @brief Updates the tree with new point locations without reallocating memory.
!> @param[inout] tree The k-d tree object to be updated.
!> @param[in] points_in A (NDims, N) array of N points with updated locations.
!> @param[in] n_active Optional: Number of active points to use (default: all points in array).
!> @details This function updates the node coordinates and rebuilds the tree structure
!>          in-place without deallocating and reallocating memory. If n_active is less
!>          than or equal to the allocated capacity, no reallocation occurs.
!>          The dimensionality must match the original tree.
subroutine kdtree_update(tree, points_in)
   type(KdTreeType), intent(INOUT) :: tree
   real(ReKi), dimension(:, :), intent(IN) :: points_in
   integer(IntKi) :: i, n_pts, n_capacity, n_dims

   n_dims = size(points_in, 1)
   n_pts = size(points_in, 2)
   
   ! Verify that the tree has been initialized
   if (.not. allocated(tree%Nodes)) then
      ! Tree is not initialized, call build instead
      call kdtree_build(tree, points_in)
      return
   end if
   
   ! Check allocated capacity and dimensions
   n_capacity = size(tree%Nodes, 2) - 1
   
   ! Check that dimensionality matches
   if (n_dims /= size(tree%Nodes, 1)) then
      ! Dimensionality changed, need to rebuild with reallocation
      call kdtree_destroy(tree)
      call kdtree_build(tree, points_in, n_max=n_capacity)
      return
   end if
   
   if (n_pts > n_capacity) then
      ! Number of points exceeds capacity, need to rebuild with reallocation
      call kdtree_destroy(tree)
      call kdtree_build(tree, points_in, n_max=n_pts)
      return
   end if
   
   ! Update number of active nodes
   tree%NNodes = n_pts

   ! Update node coordinates for active nodes
   tree%Nodes(:, 1:n_pts) = points_in(:, 1:n_pts)
   
   ! Reset original indices to sequential order for active nodes
   tree%OriginalIndices(1:n_pts) = [(i, i=1, n_pts)]
   
   ! Reset dimensions for active nodes
   tree%Dimensions(1:n_pts) = 0
   
   ! Rebuild the tree structure in-place using only active nodes
   call build_recursive(tree, 1, n_pts, 0)
end subroutine kdtree_update

!> @brief Recursively builds the tree by partitioning array segments.
recursive subroutine build_recursive(tree, start_idx, end_idx, depth)
   type(KdTreeType), intent(INOUT) :: tree
   integer, intent(IN) :: start_idx, end_idx, depth
   integer(IntKi) :: dim, mid_idx, n_dims

   if (start_idx > end_idx) return

   n_dims = size(tree%Nodes, 1)
   dim = MOD(depth, n_dims) + 1 ! Cycle through all dimensions
   mid_idx = start_idx + (end_idx - start_idx)/2

   ! Find the median of the current slice and partition around it.
   call find_median_and_partition(tree, start_idx, end_idx, mid_idx, dim)

   tree%Dimensions(mid_idx) = dim

   ! Recursively build the left and right subtrees.
   call build_recursive(tree, start_idx, mid_idx - 1, depth + 1)
   call build_recursive(tree, mid_idx + 1, end_idx, depth + 1)
end subroutine build_recursive

!> @brief Uses quickselect to find the k-th element (median) and partition the array.
subroutine find_median_and_partition(tree, left_in, right_in, k, dim)
   type(KdTreeType), intent(INOUT) :: tree
   integer, intent(IN) :: left_in, right_in, k, dim
   integer(IntKi) :: left, right
   integer(IntKi) :: p_idx

   left = left_in
   right = right_in

   do
      if (left == right) exit
      p_idx = partition(tree, left, right, dim)
      if (k == p_idx) then
         exit
      else if (k < p_idx) then
         right = p_idx - 1
      else
         left = p_idx + 1
      end if
   end do
end subroutine find_median_and_partition

!> @brief Partitions a sub-array around a pivot (last element). Lomuto partition scheme.
function partition(tree, left, right, dim) result(p_idx)
   type(KdTreeType), intent(INOUT) :: tree
   integer, intent(IN) :: left, right, dim
   integer(IntKi) :: p_idx
   real(ReKi) :: pivot_value
   integer(IntKi) :: i, store_idx

   ! For simplicity, we use the rightmost element as the pivot.
   ! A more robust implementation might use median-of-three to choose the pivot.
   pivot_value = tree%Nodes(dim, right)
   store_idx = left

   do i = left, right - 1
      if (tree%Nodes(dim, i) < pivot_value) then
         call swap_nodes(tree, i, store_idx)
         store_idx = store_idx + 1
      end if
   end do

   call swap_nodes(tree, store_idx, right)
   p_idx = store_idx
end function partition

!> @brief Swaps two nodes (and their original indices) in the tree arrays.
subroutine swap_nodes(tree, i, j)
   type(KdTreeType), intent(INOUT) :: tree
   integer, intent(IN) :: i, j
   integer(IntKi) :: temp_idx

   tree%Nodes(:, 0) = tree%Nodes(:, i)
   tree%Nodes(:, i) = tree%Nodes(:, j)
   tree%Nodes(:, j) = tree%Nodes(:, 0)

   temp_idx = tree%OriginalIndices(i)
   tree%OriginalIndices(i) = tree%OriginalIndices(j)
   tree%OriginalIndices(j) = temp_idx
end subroutine swap_nodes

!> @brief Recursively searches for the nearest neighbor.
recursive subroutine search_recursive(tree, start_idx, end_idx, query_point, best_dist_sq, best_idx)
   type(KdTreeType), intent(IN) :: tree
   integer, intent(IN) :: start_idx, end_idx
   real(ReKi), dimension(:), intent(IN) :: query_point
   real(ReKi), intent(INOUT) :: best_dist_sq
   integer, intent(INOUT) :: best_idx
   integer(IntKi) :: mid_idx, dim, i
   real(ReKi) :: d_sq, dx
   integer(IntKi) :: good_side_start, good_side_end
   integer(IntKi) :: bad_side_start, bad_side_end

   if (start_idx > end_idx) return

   ! The root of the current subtree is at its median index.
   mid_idx = start_idx + (end_idx - start_idx)/2

   ! 1. Check current node against the best so far.
   d_sq = 0.0_ReKi
   do i = 1, size(tree%Nodes, 1)
      d_sq = d_sq + (tree%Nodes(i, mid_idx) - query_point(i))**2
   end do

   if (d_sq < best_dist_sq) then
      best_dist_sq = d_sq
      best_idx = tree%OriginalIndices(mid_idx)
   end if

   dim = tree%Dimensions(mid_idx)

   ! 2. Decide which subtree to search first.
   if (query_point(dim) < tree%Nodes(dim, mid_idx)) then
      good_side_start = start_idx
      good_side_end = mid_idx - 1
      bad_side_start = mid_idx + 1
      bad_side_end = end_idx
   else
      good_side_start = mid_idx + 1
      good_side_end = end_idx
      bad_side_start = start_idx
      bad_side_end = mid_idx - 1
   end if

   ! 3. Recurse down the "good" side.
   call search_recursive(tree, good_side_start, good_side_end, query_point, best_dist_sq, best_idx)

   ! 4. Check if the "bad" side needs to be searched.
   ! This is the core of the k-d search optimization.
   dx = query_point(dim) - tree%Nodes(dim, mid_idx)
   if (dx**2 < best_dist_sq) then
      call search_recursive(tree, bad_side_start, bad_side_end, query_point, best_dist_sq, best_idx)
   end if
end subroutine search_recursive

!> @brief Recursively searches for all points within a given radius.
recursive subroutine radius_search_recursive(tree, start_idx, end_idx, query_point, radius_sq, indices, n_found)
   type(KdTreeType), intent(IN) :: tree
   integer, intent(IN) :: start_idx, end_idx
   real(ReKi), dimension(:), intent(IN) :: query_point
   real(ReKi), intent(IN) :: radius_sq
   integer, dimension(:), intent(INOUT) :: indices
   integer, intent(INOUT) :: n_found
   integer(IntKi) :: mid_idx, dim, i
   real(ReKi) :: d_sq, dx
   integer(IntKi) :: good_side_start, good_side_end
   integer(IntKi) :: bad_side_start, bad_side_end

   if (start_idx > end_idx) return

   ! The root of the current subtree is at its median index.
   mid_idx = start_idx + (end_idx - start_idx)/2

   ! 1. Check if current node is within radius.
   d_sq = 0.0_ReKi
   do i = 1, size(tree%Nodes, 1)
      d_sq = d_sq + (tree%Nodes(i, mid_idx) - query_point(i))**2
   end do

   if (d_sq <= radius_sq) then
      n_found = n_found + 1
      indices(n_found) = tree%OriginalIndices(mid_idx)
   end if

   dim = tree%Dimensions(mid_idx)

   ! 2. Decide which subtree to search first.
   if (query_point(dim) < tree%Nodes(dim, mid_idx)) then
      good_side_start = start_idx
      good_side_end = mid_idx - 1
      bad_side_start = mid_idx + 1
      bad_side_end = end_idx
   else
      good_side_start = mid_idx + 1
      good_side_end = end_idx
      bad_side_start = start_idx
      bad_side_end = mid_idx - 1
   end if

   ! 3. Recurse down the "good" side.
   call radius_search_recursive(tree, good_side_start, good_side_end, query_point, radius_sq, indices, n_found)

   ! 4. Check if the "bad" side needs to be searched.
   dx = query_point(dim) - tree%Nodes(dim, mid_idx)
   if (dx**2 <= radius_sq) then
      call radius_search_recursive(tree, bad_side_start, bad_side_end, query_point, radius_sq, indices, n_found)
   end if
end subroutine radius_search_recursive

!> @brief Sorts indices using insertion sort algorithm.
!> @details Efficient for small arrays (typically < 20 elements).
subroutine insertion_sort(indices, left, right)
   integer, dimension(:), intent(INOUT) :: indices
   integer, intent(IN) :: left, right
   integer :: i, j, key
   
   do i = left + 1, right
      key = indices(i)
      j = i - 1
      do while (j >= left)
         if (indices(j) <= key) exit
         indices(j + 1) = indices(j)
         j = j - 1
      end do
      indices(j + 1) = key
   end do
end subroutine insertion_sort

!> @brief Sorts indices using quicksort algorithm.
!> @details Recursive quicksort with insertion sort for small partitions.
recursive subroutine quicksort_indices(indices, left, right)
   integer, dimension(:), intent(INOUT) :: indices
   integer, intent(IN) :: left, right
   integer :: pivot_idx
   
   if (right - left < 20) then
      ! Use insertion sort for small partitions
      call insertion_sort(indices, left, right)
      return
   end if
   
   if (left < right) then
      pivot_idx = partition_indices(indices, left, right)
      call quicksort_indices(indices, left, pivot_idx - 1)
      call quicksort_indices(indices, pivot_idx + 1, right)
   end if
end subroutine quicksort_indices

!> @brief Partitions array for quicksort using median-of-three pivot selection.
function partition_indices(indices, left, right) result(p_idx)
   integer, dimension(:), intent(INOUT) :: indices
   integer, intent(IN) :: left, right
   integer :: p_idx
   integer :: pivot_value, i, store_idx, mid, temp
   
   ! Median-of-three pivot selection for better performance
   mid = left + (right - left) / 2
   
   ! Sort left, mid, right
   if (indices(mid) < indices(left)) then
      temp = indices(left)
      indices(left) = indices(mid)
      indices(mid) = temp
   end if
   if (indices(right) < indices(left)) then
      temp = indices(left)
      indices(left) = indices(right)
      indices(right) = temp
   end if
   if (indices(right) < indices(mid)) then
      temp = indices(mid)
      indices(mid) = indices(right)
      indices(right) = temp
   end if
   
   ! Use middle element as pivot
   pivot_value = indices(mid)
   
   ! Move pivot to end
   temp = indices(mid)
   indices(mid) = indices(right - 1)
   indices(right - 1) = temp
   
   store_idx = left
   do i = left, right - 2
      if (indices(i) < pivot_value) then
         temp = indices(i)
         indices(i) = indices(store_idx)
         indices(store_idx) = temp
         store_idx = store_idx + 1
      end if
   end do
   
   ! Move pivot to its final position
   temp = indices(store_idx)
   indices(store_idx) = indices(right - 1)
   indices(right - 1) = temp
   
   p_idx = store_idx
end function partition_indices

end module
