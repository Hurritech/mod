#ifndef MOD_PY_COMMON_HPP
#define MOD_PY_COMMON_HPP

// TODO: https://github.com/boostorg/python/pull/296
#define BOOST_BIND_GLOBAL_PLACEHOLDERS

#include <boost/python.hpp>

#undef BOOST_BIND_GLOBAL_PLACEHOLDERS

#include <optional>
#include <tuple>

namespace py = boost::python;
namespace mod::Py {

struct AttributeIsNotReadable {
};

inline void noGet(AttributeIsNotReadable) {}

// From https://wiki.python.org/moin/boost.python/HowTo#Dynamic_template_to-python_converters
template<class T1, class T2>
struct PairToTupleConverter {
	static PyObject *convert(const std::pair<T1, T2> &pair) {
		return py::incref(py::make_tuple(pair.first, pair.second).ptr());
	}
};

template<class T1, class T2, class T3>
struct TripleToTupleConverter {
	static PyObject *convert(const std::tuple<T1, T2, T3> &triple) {
		return py::incref(py::make_tuple(std::get<0>(triple), std::get<1>(triple), std::get<2>(triple)).ptr());
	}
};

// From https://stackoverflow.com/questions/16497889/how-to-expose-stdpair-to-python-using-boostpython
template<typename T1, typename T2>
struct TupleToPairConverter {
	TupleToPairConverter() {
		py::converter::registry::push_back(&convertible, &construct, py::type_id<std::pair<T1, T2>>());
	}

	static void *convertible(PyObject *obj) {
		if(!PyTuple_CheckExact(obj)) return nullptr;
		if(PyTuple_Size(obj) != 2) return nullptr;
		return obj;
	}

	static void construct(PyObject *obj, py::converter::rvalue_from_python_stage1_data *data) {
		py::tuple tuple(py::borrowed(obj));
		void *storage = reinterpret_cast<py::converter::rvalue_from_python_storage<std::pair<T1, T2>> *>(data)->storage.bytes;
		storage = new(storage) std::pair<T1, T2>(py::extract<T1>(tuple[0]), py::extract<T2>(tuple[1]));
		data->convertible = storage;
	}
};

// https://stackoverflow.com/questions/36485840/wrap-boostoptional-using-boostpython
template<typename T>
struct ToPythonOptionalValue {
	static PyObject *convert(std::optional<T> obj) {
		return obj
		       ? py::incref(py::object(*obj).ptr())
		       : py::incref(py::object().ptr());
	}
};

} // namespace mod::Py

#endif // MOD_PY_COMMON_HPP
