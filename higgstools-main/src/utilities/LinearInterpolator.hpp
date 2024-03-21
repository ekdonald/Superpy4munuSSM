/**
 * @file LinearInterpolator.hpp
 * @author Jonas Wittbrodt (jonas.wittbrodt@desy.de)
 *
 * @brief Multidimensional linear interpolation
 *
 * @copyright Copyright 2020 by the authors.
 * This file is part of HiggsBounds.
 * HiggsBounds is released under the GPLv3+.
 */
#pragma once

#include "utilities/Concepts.hpp"
#include <algorithm>
#include <array>
#include <cstddef>
#include <functional>
#include <range/v3/algorithm/any_of.hpp>
#include <range/v3/algorithm/copy.hpp>
#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/algorithm/is_sorted.hpp>
#include <range/v3/algorithm/transform.hpp>
#include <range/v3/algorithm/upper_bound.hpp>
#include <range/v3/functional/comparisons.hpp>
#include <range/v3/functional/identity.hpp>
#include <range/v3/numeric/accumulate.hpp>
#include <range/v3/numeric/inner_product.hpp>
#include <range/v3/view/cartesian_product.hpp>
#include <range/v3/view/concat.hpp>
#include <range/v3/view/iota.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/zip.hpp>
#include <range/v3/view/zip_with.hpp>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace Higgs::utilities {

namespace Detail {
/**
 * Returns an index `i` such that the closed interval `[edge[i], edge[i+1]]` is
 * valid and contains the value. Valid means that `i+1 < edge.size()`. *The
 * value has to lie in the range of the edge.*
 * @param edge the grid edge to use
 * @param value the value to place
 * @return index on the edge
 */
inline auto findAnchor(const std::vector<double> &edge, double value) noexcept {
    const auto bound = ranges::upper_bound(edge, value) - 1;
    const auto loc = static_cast<std::size_t>(bound - edge.begin());
    return bound == edge.end() - 1 ? loc - 1 : loc;
}
} // namespace Detail

/**
 * @brief A view of data in a regular grid.
 *
 * **Only stores references to the underlying data, if these are invalidated you
 * get undefined behaviour.**
 *
 * @tparam Value value type of the data
 * @tparam dim dimensionality of thr grid
 */
template <typename Value, std::size_t dim> class RegularGridView {
    static_assert(dim > 0, "A grid needs a positive dimensionality.");

  public:
    //! type of the grid edges
    using EdgeType = std::vector<double>;
    //! type for the whole grid
    using GridType = std::array<EdgeType, dim>;
    //! type of the (flattened) data
    using DataType = std::vector<Value>;
    //! type of a point in the grid
    using PointType = std::array<double, dim>;
    //!< type of an index to a point in the grid
    using IndexType = std::array<std::size_t, dim>;

    //! Access the value at the specified grid index. If the value pointer is
    //! null or the index is out of bounds return `Value{}`.
    const Value &operator[](const IndexType &index) const noexcept {
        std::size_t flatIndex = 0;
        for (std::size_t k = 0; k != dim; ++k) {
            std::size_t product = index[k];
            for (std::size_t l = k + 1; l != dim; ++l) {
                product *= sizes_[l];
            }
            flatIndex += product;
        }
        if (values_ == nullptr || flatIndex >= values_->size()) {
            static auto defVal = Value{}; // LCOV_EXCL_LINE
            return defVal;
        }
        return (*values_)[flatIndex];
    }

    //! A point on the grid expressed through the index of an anchor grid point
    //! and normalizes coordinates in the unit cell starting at the anchor.
    struct LocatedPoint {
        //! Index of the anchor point in the grid.
        IndexType anchor;
        //! Normalized coordinates of the point in the unit cell. In these
        //! coordinates the anchor is at all 0.
        PointType normDist;
    };

    //! Locate the point in the grid and return it as a LocatedPoint. If the
    //! point is outside the grid it is clamped to the grid extent. Returns
    //! `LocatedPoint{}` if the edge pointer is null.
    LocatedPoint locate(const PointType &point) const noexcept {
        auto result = LocatedPoint{};
        if (edges_ != nullptr) {
            auto clampedPoint = clamp(point);
            constexpr auto locateAlong = [](const EdgeType &edge, double v) {
                if (edge.size() < 2) {
                    return std::pair{0ul, 0.};
                }
                auto idx = Detail::findAnchor(edge, v);
                return std::pair{idx,
                                 (v - edge[idx]) / (edge[idx + 1] - edge[idx])};
            };
            ranges::transform(
                *edges_, clampedPoint,
                ranges::views::zip(result.anchor, result.normDist).begin(),
                locateAlong);
        }
        return result;
    }

    //! (min, max) of the grid for every dim
    const std::array<std::pair<double, double>, dim> &extent() const noexcept {
        return extent_;
    }

    //! clamp the point to the grid extent
    PointType clamp(PointType point) const noexcept {
        ranges::transform(
            point, extent_, point.begin(), [](auto val, const auto &bounds) {
                return std::clamp(val, bounds.first, bounds.second);
            });
        return point;
    }

    /**
     * @brief Construct a new Regular Grid View object.
     *
     * Only stores references to the passes edges and values. There has to be
     * exactly one value for each point in the grid spanned by the edges and all
     * edges must be sorted.
     *
     * @param edges grid edges
     * @param values row-major order flattened grid values
     */
    RegularGridView(const GridType &edges, const DataType &values)
        : edges_{&edges}, values_{&values}, sizes_{sizes(edges)},
          extent_{extent(edges)} {
        if (ranges::accumulate(sizes_, std::size_t{1}, std::multiplies{}) !=
            values_->size()) {
            throw std::out_of_range("Size of values does not match the grid.");
        }
        ranges::for_each(*edges_, [](const EdgeType &e) {
            if (!ranges::is_sorted(e)) {
                throw std::runtime_error("Grid edges have to be sorted.");
            }
        });
    }

    //! Default constructor.
    RegularGridView() = default;

  private:
    //! get the size of every grid dimension
    static auto sizes(const GridType &edges) noexcept {
        auto result = std::array<size_t, dim>{};
        ranges::transform(edges, result.begin(),
                          [](const auto &vec) { return vec.size(); });
        return result;
    }
    //! get (min,max) of every grid dimension
    static auto extent(const GridType &edges) noexcept {
        constexpr auto getFrontBack = [](const auto &edge) {
            return std::pair{edge.front(), edge.back()};
        }; // LCOV_EXCL_LINE
        auto extent = std::array<std::pair<double, double>, dim>{};
        ranges::transform(edges, extent.begin(), getFrontBack);
        return extent;
    }

    const GridType *edges_;
    const DataType *values_;
    std::array<std::size_t, dim> sizes_ = {};
    std::array<std::pair<double, double>, dim> extent_ = {};
};

/**
 * @brief Data on a regular grid.
 *
 * Same functionality as RegularGridView, but stores its data internally and is
 * therefore much safer to use.
 *
 * @tparam Value the value type of the data
 * @tparam dim dimensionality of thr grid
 */
template <typename Value, std::size_t dim>
class RegularGrid : public RegularGridView<Value, dim> {
  public:
    //! type for the whole grid
    using GridType =
        std::array<typename RegularGridView<Value, dim>::EdgeType, dim>;
    //! type of the (flattened) data
    using DataType = std::vector<Value>;
    /**
     * @brief Construct a new Regular Grid object.
     * Copies the passed edges and values into internal storage.
     * See also RegularGridView::RegularGridView.
     * @param edges grid edges
     * @param values row-major order flattened grid values
     */
    RegularGrid(const GridType &edges, const DataType &values)
        : RegularGridView<Value, dim>{}, edges_{edges}, values_{values} {
        static_cast<RegularGridView<Value, dim> &>(*this) =
            RegularGridView<Value, dim>(edges_, values_);
    }
    //! default constructs an empty grid
    RegularGrid() = default;

    //! Copy constructor
    RegularGrid(const RegularGrid<Value, dim> &other)
        : RegularGrid{other.edges_, other.values_} {}

    //! Move constructor
    RegularGrid(RegularGrid<Value, dim> &&other)
        : RegularGridView<Value, dim>{}, edges_{std::move(other.edges_)},
          values_{std::move(other.values_)} {
        static_cast<RegularGridView<Value, dim> &>(*this) =
            RegularGridView<Value, dim>(edges_, values_);
    }

    //! copy assignment
    RegularGrid &operator=(const RegularGrid<Value, dim> &other) {
        if (this != &other) {
            edges_ = other.edges_;
            values_ = other.values_;
            static_cast<RegularGridView<Value, dim> &>(*this) =
                RegularGridView<Value, dim>(edges_, values_);
        }
        return *this;
    }

    //! move assignment
    RegularGrid &operator=(RegularGrid<Value, dim> &&other) {
        if (this != &other) {
            edges_ = std::move(other.edges_);
            values_ = std::move(other.values_);
            static_cast<RegularGridView<Value, dim> &>(*this) =
                RegularGridView<Value, dim>(edges_, values_);
        }
        return *this;
    }

    ~RegularGrid() = default;

  private:
    GridType edges_;
    DataType values_;
};

namespace Detail {

//! 2 to the j
constexpr size_t pow2(size_t j) { return 1ul << j; }

/**
 * Generates the corners of the dim-dimensional unit cell. The resulting
 * array is lexicographically sorted.
 * @tparam dim dimensionality
 * @return constexpr array of unit hypercube corners with `2^dim` elements,
 * each of size dim.
 */
template <size_t dim> constexpr auto unitCellCorners() {
    auto result = std::array<std::array<std::size_t, dim>, pow2(dim)>{};
    for (size_t i = 0; i != pow2(dim); ++i) {
        for (size_t j = 0; j != dim; ++j) {
            result[i][dim - 1 - j] = (i / pow2(j)) % 2;
        }
    }
    return result;
}

//! Equals T::dimensionality if present, zero otherwise. Default case zero.
template <typename T, typename = const std::size_t>
constexpr std::size_t valueDim = 0;
//! Equals T::dimensionality if present, zero otherwise. Specialization if
//! present.
template <typename T>
constexpr std::size_t valueDim<T, decltype(T::dimensionality)> =
    T::dimensionality;
} // namespace Detail

/**
 * Multilinear interpolation on a regular grid.
 *
 * This class takes a regular grid defined by its `dim` edges and the
 * corresponding data to construct a RegularGrid (or RegularGridView, see
 * below). Multilinear interpolation to any point within the range of the edges
 * can then be performed using operator()(). Performs nearest neighbor
 * extrapolation if a requested value lies outside the grid (i.e. clamps the
 * value to the grid edges.)
 *
 * Since multilinear interpolation scales with `2**dim` this is not the most
 * efficient method for very large `dim`. The algorithm used is the same as in
 * `scipy.interpolate.RegularGridInterpolator`. Interpolation with
 * multidimensional values is supported, by passing a vector of multidimensional
 * objects as values (see ValueType).
 *
 * It is also possible to use another LinearInterpolator as ValueType. In this
 * case, any point passed to operator()() is split, with the first dim entries
 * used for the current interpolation and the remainder passed to the values at
 * the grid edges.
 *
 * @tparam dim the dimensionality of the underlying regular grid
 * @tparam ValueType the std::default_initializable type of the values on the
 * grid
 * @tparam view use a RegularGridView to access the data instead of copying it
 * into a RegularGrid
 */
template <size_t dim, class Value = double, bool view = false>
REQUIRES(std::default_initializable<Value>)
class LinearInterpolator {
    //! type of the underlying RegularGrid (it's a RegularGridView if
    //! `view==true`)
    using RegularGridType =
        std::conditional_t<view, RegularGridView<Value, dim>,
                           RegularGrid<Value, dim>>;

  public:
    //! the complete dimensionality of this interpolator, including any
    //! dimensionality the Value%s themselves might have
    static constexpr std::size_t dimensionality = dim + Detail::valueDim<Value>;
    //! type for the whole grid of this interpolator
    using GridType = std::array<typename RegularGridType::EdgeType, dim>;
    //! type of the (flattened) data
    using DataType = std::vector<Value>;
    //! type of a point in the complete grid
    using PointType = std::array<double, dimensionality>;

    //! Return the extent of the underlying grid.
    const std::array<std::pair<double, double>, dim> &extent() const noexcept {
        return grid_.extent();
    }

    //! Obtain a multilinear interpolated value at the given point.
    auto operator()(const PointType &point) const noexcept {
        if constexpr (dim == dimensionality) {
            const auto contribution =
                [this, loc = grid_.locate(point)](const auto &corner) {
                    return grid_[cellIndex(loc.anchor, corner)] *
                           cornerWeight(corner, loc.normDist);
                };
            return ranges::accumulate(corners_, Value{}, std::plus{},
                                      contribution);
        } else {
            auto pointParts = splitPoint(point);
            const auto contribution = [this,
                                       loc = grid_.locate(pointParts.first),
                                       &pointParts](const auto &corner) {
                return grid_[cellIndex(loc.anchor, corner)](pointParts.second) *
                       cornerWeight(corner, loc.normDist);
            };
            return ranges::accumulate(
                corners_,
                std::invoke_result_t<
                    Value, std::array<double, dimensionality - dim>>{},
                std::plus{}, contribution);
        }
    }

    //! Returns the maximum value in the hypercube spanned by
    //! the two given corners. Includes the values at all corners of the cube as
    //! well as at all interior grid points.
    template <class Comparator = std::less<>, class Projector = ranges::identity>
    Value maxWithin(const PointType &corner1, const PointType &corner2,
                    const Comparator &comparator = Comparator{},
                    const Projector &projector = Projector{}) const noexcept {
        return ranges::max(
            ranges::views::concat(atAllCorners(corner1, corner2),
                                  atAllPointsBetween(corner1, corner2)),
            comparator, projector);
    }

    //! Returns the minimum value in the hypercube spanned by
    //! the two given corners. Includes the values at all corners of the cube as
    //! well as at all interior grid points.
    template <class Comparator = std::less<>, class Projector = ranges::identity>
    Value minWithin(const PointType &corner1, const PointType &corner2,
                    const Comparator &comparator = Comparator{},
                    const Projector &projector = Projector{}) const noexcept {
        return ranges::min(
            ranges::views::concat(atAllCorners(corner1, corner2),
                                  atAllPointsBetween(corner1, corner2)),
            comparator, projector);
    }

    //! Construct a new Linear Interpolator object. Passes the given
    //! interpolation grid parameters to the underlying RegularGrid (or
    //! RegularGridView, if `view` is true).
    LinearInterpolator(const GridType &edges, const DataType &values)
        : grid_{edges, values} {}

    //! Default constructs a LinearInterpolator. Returns Value{} for any call to
    //! operator()().
    LinearInterpolator() = default;

    //! Return the multilinear interpolated values at all corners of the box
    //! spanned by the two given opposing corners.
    auto atAllCorners(const PointType &corner1,
                      const PointType &corner2) const noexcept {
        static constexpr auto corners =
            Detail::unitCellCorners<dimensionality>();
        const auto fromUnitCorner = [&corner1, &corner2](const auto &uc) {
            constexpr auto selectElement = [](std::size_t uce, double e1,
                                              double e2) {
                return uce == 0 ? e1 : e2;
            };
            auto res = PointType{};
            ranges::copy(
                ranges::views::zip_with(selectElement, uc, corner1, corner2),
                res.begin());
            return res;
        };
        return corners | ranges::views::transform(fromUnitCorner) |
               ranges::views::transform(*this);
    }

    //! Returns the values at all grid points fully enclosed inside the
    //! hypercube spanned by the two given corners.
    auto atAllPointsBetween(const PointType &corner1,
                            const PointType &corner2) const {
        static_assert(dim == dimensionality);
        auto locatedCorner1 = grid_.locate(corner1);
        auto locatedCorner2 = grid_.locate(corner2);
        return allPointsBetweenImpl(locatedCorner1.anchor,
                                    locatedCorner2.anchor,
                                    std::make_index_sequence<dim>{}) |
               ranges::views::transform([this](const auto &tuple) {
                   constexpr auto asArray = [](auto &&...x) {
                       return std::array{x...};
                   };
                   return grid_[std::apply(asArray, tuple)];
               });
        ;
    }

  private:
    template <std::size_t pointDim, std::size_t... Is>
    auto
    allPointsBetweenImpl(const std::array<std::size_t, pointDim> &firstCorner,
                         const std::array<std::size_t, pointDim> &secondCorner,
                         std::index_sequence<Is...>) const {
        return ranges::views::cartesian_product(indicesBetween(
            std::get<Is>(firstCorner), std::get<Is>(secondCorner))...);
    }

    auto indicesBetween(std::size_t a, std::size_t b) const {
        return ranges::views::iota(std::min(a, b) + 1, std::max(a, b) + 1);
    }

    //! compute the weight of the given corner of the init cell
    double cornerWeight(
        const typename RegularGridType::IndexType &corner,
        const typename RegularGridType::PointType &normDist) const noexcept {
        return ranges::inner_product(
            corner, normDist, 1., std::multiplies{},
            [](auto c, auto d) { return c == 0 ? 1 - d : d; });
    }

    //! Compute the index of the corner of a cell corner by adding a unit cell
    //! corner to an anchor.
    auto cellIndex(
        const typename RegularGridType::IndexType &anchor,
        const typename RegularGridType::IndexType &corner) const noexcept {
        auto index = typename RegularGridType::IndexType{};
        ranges::transform(
            anchor, corner, index.begin(),
            [](std::size_t i1, std::size_t i2) { return i1 + i2; });
        return index;
    }

    //! Splits a point into the part that corresponds to this interpolators grid
    //! and a remainder that is passed on the underlying value.
    auto splitPoint(const PointType &point) const noexcept {
        auto first = typename RegularGridType::PointType{};
        std::copy(point.begin(), point.begin() + dim, first.begin());
        auto second = std::array<double, dimensionality - dim>{};
        std::copy(point.begin() + dim, point.end(), second.begin());
        return std::pair{first, second};
    }

    static constexpr auto corners_ = Detail::unitCellCorners<dim>();
    RegularGridType grid_;
};

} // namespace Higgs::utilities
