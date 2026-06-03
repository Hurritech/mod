#ifndef MOD_CAUSALITY_STOCHSIM_HPP
#define MOD_CAUSALITY_STOCHSIM_HPP

#include <mod/BuildConfig.hpp>
#include <mod/causality/ForwardDecl.hpp>
#include <mod/causality/EventTrace.hpp>
#include <mod/dg/DG.hpp>
#include <mod/dg/GraphInterface.hpp>

#include <optional>
#include <tuple>

namespace mod::causality {

// rst-class: causality::DrawMassActionFunction
// rst:
// rst:		A helper class for performing stochastic simulations where
// rst:		events are drawn according to the law of mass action.
// rst:		Importantly, if the underlying derivation graph is enlarged then
// rst:		:func:`syncSize` must be called before calling :func:`draw`.
// rst:
// rst:		.. versionadded:: 1.1
// rst:
// rst:		.. todo:: check version added before release
// rst:
// rst-class-start:
struct MOD_DECL DrawMassActionFunction {
	// rst: .. function:: DrawMassActionFunction(std::shared_ptr<dg::DG> dg_,\
	// rst:                                      std::function<std::pair<double, bool>(dg::DG::Vertex)> inputRate, \
	// rst:                                      std::function<std::pair<double, bool>(dg::DG::HyperEdge)> reactionRate, \
	// rst:                                      std::function<std::pair<double, bool>(dg::DG::Vertex)> outputRate)
	// rst:
	// rst:		Construct a new instance, based on the given derivation graph.
	// rst:		The return value of the rate callbacks must be 1) the rate and 2) a boolean telling whether
	// rst:		the library should cache the rate. If `true` no more calls with the same argument will be made.
	// rst:		Each of the rate functions may be an empty `std::function` (i.e., default constructed),
	// rst:		which means a default rate is used: input rate 0.0, reaction rate 1.0, output rate 0.0.
	// rst:
	// rst:		:throws: :class:`LogicError` if `!dg_`.
	// rst:		:throws: :class:`LogicError` if neither `dg_->hasActiveBuilder()` nor `dg_->isLocked()`.
	DrawMassActionFunction(std::shared_ptr<dg::DG> dg_,
	                   std::function<std::pair<double, bool>(dg::DG::Vertex)> inputRate,
	                   std::function<std::pair<double, bool>(dg::DG::HyperEdge)> reactionRate,
	                   std::function<std::pair<double, bool>(dg::DG::Vertex)> outputRate);
	~DrawMassActionFunction();
	DrawMassActionFunction(DrawMassActionFunction &&);
	DrawMassActionFunction &operator=(DrawMassActionFunction &&);
	DrawMassActionFunction(const DrawMassActionFunction &);
	DrawMassActionFunction &operator=(const DrawMassActionFunction &);
	// rst: .. function:: void syncSize()
	// rst:
	// rst:		Enlarges the internal data structures to the current size of the underlying derivation graph.
	void syncSize();
	// rst: .. function:: std::tuple<std::optional<Action>, double, bool> draw(const Marking &m)
	// rst:
	// rst:		:returns: an action, randomly selected according to the low of mass action
	// rst:			based on the given population/marking and the rates.
	// rst:			The second returned component is the sum of reactivity for the given reactions/hyperedges.
	// rst:         The third returned component is false to indicate that no time increment is returned
	// rst:
	// rst:			If no events are possible, then the second component is 0.0 and the first component has no value.
	// rst:
	// rst:			Use :func:`rngReseed` to seed the pseudo-random bit generator used for the selection.
	// rst:
	// rst:		:throws: :class:`LogicError` if `m` is not a marking on the underlying derivation graph.
	// rst:
	// rst:		Requires :func:`syncSize` to have been called since the last time the underlying derivation graph has changed size.
	std::tuple<std::optional<Action>, double, bool> draw(const Marking &m);
private:
	struct Pimpl;
	std::unique_ptr<Pimpl> p;
};
// rst-class-end:

// ==============================================================================================

// rst-class: causality::DrawMassActionTauLeapingFunction
// rst:
// rst:		A helper class for performing stochastic simulations where
// rst:		events are drawn according to the law of mass action, but
// rst:     batched according to the Tau-Leaping technique.
// rst:		Importantly, if the underlying derivation graph is enlarged then
// rst:		:func:`syncSize` must be called before calling :func:`draw`.
// rst:
// rst:		.. versionadded:: 1.1
// rst:
// rst:		.. todo:: check version added before release
// rst:
// rst-class-start:
struct MOD_DECL DrawMassActionTauLeapingFunction {
	// rst: .. function:: DrawMassActionTauLeapingFunction(std::shared_ptr<dg::DG> dg_,\
	// rst:        std::function<std::pair<double, bool>(dg::DG::Vertex)> inputRate, \
	// rst:        std::function<std::pair<double, bool>(dg::DG::HyperEdge)> reactionRate, \
	// rst:        std::function<std::pair<double, bool>(dg::DG::Vertex)> outputRate, \
	// rst:        int dc, \
	// rst:        double epsilon)
	// rst:
	// rst:		Construct a new instance, based on the given derivation graph.
	// rst:		The return value of the rate callbacks must be 1) the rate and 2) a boolean telling whether
	// rst:		the library should cache the rate. If `true` no more calls with the same argument will be made.
	// rst:		Each of the rate functions may be an empty `std::function` (i.e., default constructed),
	// rst:		which means a default rate is used: input rate 0.0, reaction rate 1.0, output rate 0.0.
	// rst:
	// rst:		:throws: :class:`LogicError` if `!dg_`.
	// rst:		:throws: :class:`LogicError` if neither `dg_->hasActiveBuilder()` nor `dg_->isLocked()`.
	DrawMassActionTauLeapingFunction(std::shared_ptr<dg::DG> dg_,
	                   std::function<std::pair<double, bool>(dg::DG::Vertex)> inputRate,
	                   std::function<std::pair<double, bool>(dg::DG::HyperEdge)> reactionRate,
	                   std::function<std::pair<double, bool>(dg::DG::Vertex)> outputRate,
	                   int dc,
	                   double epsilon);
	~DrawMassActionTauLeapingFunction();
	DrawMassActionTauLeapingFunction(DrawMassActionTauLeapingFunction &&);
	DrawMassActionTauLeapingFunction &operator=(DrawMassActionTauLeapingFunction &&);
	DrawMassActionTauLeapingFunction(const DrawMassActionTauLeapingFunction &);
	DrawMassActionTauLeapingFunction &operator=(const DrawMassActionTauLeapingFunction &);
	// rst: .. function:: void syncSize()
	// rst:
	// rst:		Enlarges the internal data structures to the current size of the underlying derivation graph.
	void syncSize();
	// rst: .. function:: std::tuple<std::optional<Action>, double, bool> draw(const Marking &m)
	// rst:
	// rst:		:returns: an action, randomly selected according to the law of mass action
	// rst:			based on the given population/marking and the rates.
	// rst:			The second returned component is the tau (time increment) with the current update.
	// rst:         The third returned component is true to indicate that a time increment is returned
	// rst:
	// rst:			If no events are possible, then the second component is 0.0 and the first component has no value.
	// rst:
	// rst:			Use :func:`rngReseed` to seed the pseudo-random bit generator used for the selection.
	// rst:
	// rst:		:throws: :class:`LogicError` if `m` is not a marking on the underlying derivation graph.
	// rst:
	// rst:		Requires :func:`syncSize` to have been called since the last time the underlying derivation graph has changed size.
	std::tuple<std::optional<Action>, double, bool> draw(const Marking &m);
private:
	struct Pimpl;
	std::unique_ptr<Pimpl> p;
};
// rst-class-end:

// ==============================================================================================

// rst-class: causality::DrawMassActionEulerMaruyamaFunction
// rst:
// rst:		A helper class for performing stochastic simulations where
// rst:		the state is updates according to the Euler-Maruyama-method for
// rst:     stochastic differential equations.
// rst:		Importantly, if the underlying derivation graph is enlarged then
// rst:		:func:`syncSize` must be called before calling :func:`draw`.
// rst:
// rst:		.. versionadded:: 1.1
// rst:
// rst:		.. todo:: check version added before release
// rst:
// rst-class-start:
struct MOD_DECL DrawMassActionEulerMaruyamaFunction {
	// rst: .. function:: DrawMassActionEulerMaruyamaFunction(std::shared_ptr<dg::DG> dg_,\
	// rst:        std::function<std::pair<double, bool>(dg::DG::Vertex)> inputRate, \
	// rst:        std::function<std::pair<double, bool>(dg::DG::HyperEdge)> reactionRate, \
	// rst:        std::function<std::pair<double, bool>(dg::DG::Vertex)> outputRate, \
	// rst:        double tau)
	// rst:
	// rst:		Construct a new instance, based on the given derivation graph.
	// rst:		The return value of the rate callbacks must be 1) the rate and 2) a boolean telling whether
	// rst:		the library should cache the rate. If `true` no more calls with the same argument will be made.
	// rst:		Each of the rate functions may be an empty `std::function` (i.e., default constructed),
	// rst:		which means a default rate is used: input rate 0.0, reaction rate 1.0, output rate 0.0.
	// rst:
	// rst:		:throws: :class:`LogicError` if `!dg_`.
	// rst:		:throws: :class:`LogicError` if neither `dg_->hasActiveBuilder()` nor `dg_->isLocked()`.
	DrawMassActionEulerMaruyamaFunction(std::shared_ptr<dg::DG> dg_,
	                   std::function<std::pair<double, bool>(dg::DG::Vertex)> inputRate,
	                   std::function<std::pair<double, bool>(dg::DG::HyperEdge)> reactionRate,
	                   std::function<std::pair<double, bool>(dg::DG::Vertex)> outputRate,
	                   double tau);
	~DrawMassActionEulerMaruyamaFunction();
	DrawMassActionEulerMaruyamaFunction(DrawMassActionEulerMaruyamaFunction &&);
	DrawMassActionEulerMaruyamaFunction &operator=(DrawMassActionEulerMaruyamaFunction &&);
	DrawMassActionEulerMaruyamaFunction(const DrawMassActionEulerMaruyamaFunction &);
	DrawMassActionEulerMaruyamaFunction &operator=(const DrawMassActionEulerMaruyamaFunction &);
	// rst: .. function:: void syncSize()
	// rst:
	// rst:		Enlarges the internal data structures to the current size of the underlying derivation graph.
	void syncSize();
	// rst: .. function:: std::tuple<std::optional<Action>, double, bool> draw(const Marking &m)
	// rst:
	// rst:		:returns: an action, randomly selected according to the law of mass action
	// rst:			based on the given population/marking and the rates.
	// rst:			The second returned component is the tau (time increment) with the current update.
	// rst:         The third returned component is true to indicate that a time increment is returned
	// rst:
	// rst:			If no events are possible, then the second component is 0.0 and the first component has no value.
	// rst:
	// rst:			Use :func:`rngReseed` to seed the pseudo-random bit generator used for the selection.
	// rst:
	// rst:		:throws: :class:`LogicError` if `m` is not a marking on the underlying derivation graph.
	// rst:
	// rst:		Requires :func:`syncSize` to have been called since the last time the underlying derivation graph has changed size.
	std::tuple<std::optional<Action>, double, bool> draw(const Marking &m);
private:
	struct Pimpl;
	std::unique_ptr<Pimpl> p;
};
// rst-class-end:

// ==============================================================================================

struct MOD_DECL DrawMassActionSKRockFunction {
	DrawMassActionSKRockFunction(std::shared_ptr<dg::DG> dg_,
	                   std::function<std::pair<double, bool>(dg::DG::Vertex)> inputRate,
	                   std::function<std::pair<double, bool>(dg::DG::HyperEdge)> reactionRate,
	                   std::function<std::pair<double, bool>(dg::DG::Vertex)> outputRate,
	                   double tau,
	                   int stages);
	~DrawMassActionSKRockFunction();
	DrawMassActionSKRockFunction(DrawMassActionSKRockFunction &&);
	DrawMassActionSKRockFunction &operator=(DrawMassActionSKRockFunction &&);
	DrawMassActionSKRockFunction(const DrawMassActionSKRockFunction &);
	DrawMassActionSKRockFunction &operator=(const DrawMassActionSKRockFunction &);
	void syncSize();
	std::tuple<std::optional<Action>, double, bool> draw(const Marking &m);
private:
	struct Pimpl;
	std::unique_ptr<Pimpl> p;
};

// ==============================================================================================

struct MOD_DECL DrawMassActionComplexEulerMaruyamaFunction {
	DrawMassActionComplexEulerMaruyamaFunction(std::shared_ptr<dg::DG> dg_,
	                   std::function<std::pair<double, bool>(dg::DG::Vertex)> inputRate,
	                   std::function<std::pair<double, bool>(dg::DG::HyperEdge)> reactionRate,
	                   std::function<std::pair<double, bool>(dg::DG::Vertex)> outputRate,
	                   double tau);
	~DrawMassActionComplexEulerMaruyamaFunction();
	DrawMassActionComplexEulerMaruyamaFunction(DrawMassActionComplexEulerMaruyamaFunction &&);
	DrawMassActionComplexEulerMaruyamaFunction &operator=(DrawMassActionComplexEulerMaruyamaFunction &&);
	DrawMassActionComplexEulerMaruyamaFunction(const DrawMassActionComplexEulerMaruyamaFunction &);
	DrawMassActionComplexEulerMaruyamaFunction &operator=(const DrawMassActionComplexEulerMaruyamaFunction &);
	void syncSize();
	std::tuple<std::optional<Action>, double, bool> draw(const Marking &m);
private:
	struct Pimpl;
	std::unique_ptr<Pimpl> p;
};

// ==============================================================================================

struct MOD_DECL DrawMassActionComplexSKRockFunction {
	DrawMassActionComplexSKRockFunction(std::shared_ptr<dg::DG> dg_,
	                   std::function<std::pair<double, bool>(dg::DG::Vertex)> inputRate,
	                   std::function<std::pair<double, bool>(dg::DG::HyperEdge)> reactionRate,
	                   std::function<std::pair<double, bool>(dg::DG::Vertex)> outputRate,
	                   double tau,
	                   int stages);
	~DrawMassActionComplexSKRockFunction();
	DrawMassActionComplexSKRockFunction(DrawMassActionComplexSKRockFunction &&);
	DrawMassActionComplexSKRockFunction &operator=(DrawMassActionComplexSKRockFunction &&);
	DrawMassActionComplexSKRockFunction(const DrawMassActionComplexSKRockFunction &);
	DrawMassActionComplexSKRockFunction &operator=(const DrawMassActionComplexSKRockFunction &);
	void syncSize();
	std::tuple<std::optional<Action>, double, bool> draw(const Marking &m);
private:
	struct Pimpl;
	std::unique_ptr<Pimpl> p;
};

// ==============================================================================================

struct MOD_DECL SimulatorImpl {
	SimulatorImpl();
	~SimulatorImpl();
public:
	int getIteration() const;
	double getTime() const;
	void setTime_delete(double value);
public:
	void doIteration();
private:
	struct Pimpl;
	std::unique_ptr<Pimpl> p;
};

} // namespace mod::causality

#endif // MOD_CAUSALITY_STOCHSIM_HPP
