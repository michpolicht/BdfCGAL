#ifndef NASTRAN_EXPORT_HPP
#define NASTRAN_EXPORT_HPP

#include "functions.hpp"

#include <CGAL/IO/File_medit.h>

namespace CGAL {
namespace Mesh_3 {

template <typename C3T3, bool rebind, bool no_patch>
struct Nastran_pmap_generator {};


template <typename C3T3>
struct Nastran_pmap_generator<C3T3, true, false>
{
	typedef Rebind_cell_pmap<C3T3>                            Cell_pmap;
	typedef Rebind_facet_pmap<C3T3, Cell_pmap>                Facet_pmap;
	typedef Null_facet_pmap<C3T3, Cell_pmap>                  Facet_pmap_twice;
	typedef Default_vertex_pmap<C3T3, Cell_pmap, Facet_pmap>  Vertex_pmap;

	bool print_twice() {
		return false;
	}
};


template <typename C3T3>
struct Nastran_pmap_generator<C3T3, true, true>
{
	typedef Rebind_cell_pmap<C3T3>                            Cell_pmap;
	typedef No_patch_facet_pmap_first<C3T3, Cell_pmap>         Facet_pmap;
	typedef No_patch_facet_pmap_second<C3T3, Cell_pmap>        Facet_pmap_twice;
	typedef Default_vertex_pmap<C3T3, Cell_pmap, Facet_pmap>  Vertex_pmap;

	bool print_twice() {
		return true;
	}
};


template <typename C3T3>
struct Nastran_pmap_generator<C3T3, false, true>
{
	typedef No_rebind_cell_pmap<C3T3>                         Cell_pmap;
	typedef No_patch_facet_pmap_first<C3T3, Cell_pmap>         Facet_pmap;
	typedef No_patch_facet_pmap_second<C3T3, Cell_pmap>        Facet_pmap_twice;
	typedef Default_vertex_pmap<C3T3, Cell_pmap, Facet_pmap>  Vertex_pmap;

	bool print_twice() {
		return true;
	}
};

template <typename C3T3>
struct Nastran_pmap_generator<C3T3, false, false>
{
	typedef No_rebind_cell_pmap<C3T3>                         Cell_pmap;
	typedef Rebind_facet_pmap<C3T3, Cell_pmap>                 Facet_pmap;
	typedef Null_facet_pmap<C3T3, Cell_pmap>                  Facet_pmap_twice;
	typedef Null_vertex_pmap<C3T3, Cell_pmap, Facet_pmap>     Vertex_pmap;

	bool print_twice() {
		return false;
	}
};


template <class C3T3, bool rebind, bool no_patch>
void
output_to_nastran(std::ostream & os,
		const C3T3 & c3t3)
{
#ifdef CGAL_MESH_3_IO_VERBOSE
	std::cerr << "Output to nastran:\n";
#endif

	typedef Nastran_pmap_generator<C3T3, rebind, no_patch> Generator;
	typedef typename Generator::Cell_pmap Cell_pmap;
	typedef typename Generator::Facet_pmap Facet_pmap;
	typedef typename Generator::Facet_pmap_twice Facet_pmap_twice;
	typedef typename Generator::Vertex_pmap Vertex_pmap;

	Cell_pmap cell_pmap(c3t3);
	Facet_pmap facet_pmap(c3t3, cell_pmap);
	Facet_pmap_twice facet_pmap_twice(c3t3, cell_pmap);
	Vertex_pmap vertex_pmap(c3t3, cell_pmap, facet_pmap);

	output_to_nastran(os,
			c3t3,
			vertex_pmap,
			facet_pmap,
			cell_pmap,
			facet_pmap_twice,
			Generator().print_twice());

#ifdef CGAL_MESH_3_IO_VERBOSE
	std::cerr << "done.\n";
#endif
}

template <class C3T3, bool rebind, bool no_patch>
void
squashed_output_to_nastran(std::ostream & os,
		const C3T3 & c3t3,
		const std::vector<std::pair<int, int>> & squashRanges)
{
#ifdef CGAL_MESH_3_IO_VERBOSE
	std::cerr << "Output to nastran:\n";
#endif

	typedef Nastran_pmap_generator<C3T3, rebind, no_patch> Generator;
	typedef typename Generator::Cell_pmap Cell_pmap;
	typedef typename Generator::Facet_pmap Facet_pmap;
	typedef typename Generator::Facet_pmap_twice Facet_pmap_twice;
	typedef typename Generator::Vertex_pmap Vertex_pmap;

	Cell_pmap cell_pmap(c3t3);
	Facet_pmap facet_pmap(c3t3, cell_pmap);
	Facet_pmap_twice facet_pmap_twice(c3t3, cell_pmap);
	Vertex_pmap vertex_pmap(c3t3, cell_pmap, facet_pmap);

	squashed_output_to_nastran(os,
			c3t3,
			squashRanges,
			vertex_pmap,
			facet_pmap,
			cell_pmap,
			facet_pmap_twice,
			Generator().print_twice());

#ifdef CGAL_MESH_3_IO_VERBOSE
	std::cerr << "done.\n";
#endif
}


template <class C3T3,
		class Vertex_index_property_map,
		class Facet_index_property_map,
		class Facet_index_property_map_twice,
		class Cell_index_property_map>
void squashed_output_to_nastran(std::ostream & os,
		const C3T3 & c3t3,
		const std::vector<std::pair<int, int>> & squashRanges,
		const Vertex_index_property_map & vertex_pmap,
		const Facet_index_property_map & facet_pmap,
		const Cell_index_property_map & cell_pmap,
		const Facet_index_property_map_twice & facet_twice_pmap = Facet_index_property_map_twice(),
		const bool print_each_facet_twice = false)
{
	typedef typename C3T3::Triangulation Tr;
	typedef typename C3T3::Facets_in_complex_iterator Facet_iterator;
	typedef typename C3T3::Cells_in_complex_iterator Cell_iterator;

	typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
	typedef typename Tr::Vertex_handle Vertex_handle;
	typedef typename Tr::Weighted_point Weighted_point;

	const Tr & tr = c3t3.triangulation();

	//-------------------------------------------------------
	// File output
	//-------------------------------------------------------

	//-------------------------------------------------------
	// Header
	//-------------------------------------------------------
	os << std::setprecision(17);

	os << "$ Created by CGAL with MP exporter\n";


	//-------------------------------------------------------
	// Vertices
	//-------------------------------------------------------

	os << "$ Vertices: " << tr.number_of_vertices() << '\n';

	boost::unordered_map<Vertex_handle, int> V;
	int inum = 1;
	for ( Finite_vertices_iterator vit = tr.finite_vertices_begin();
			vit != tr.finite_vertices_end();
			++vit)
	{
		V[vit] = inum++;
		Weighted_point p = tr.point(vit);
		os << "GRID"
				<< ' ' << V[vit]
				<< ' ' << get(vertex_pmap, vit)
				<< ' ' << CGAL::to_double(p.x())
				<< ' ' << CGAL::to_double(p.y())
				<< ' ' << CGAL::to_double(p.z())
				<< '\n';
	}

	//-------------------------------------------------------
	// Tetrahedra
	//-------------------------------------------------------
	os << "$ Tetrahedra: " << c3t3.number_of_cells_in_complex() << '\n';

	int ctetraIndex = 1;
	for ( Cell_iterator cit = c3t3.cells_in_complex_begin() ;
			cit != c3t3.cells_in_complex_end() ;
			++cit )
	{
		os << "CTETRA";
		os << ' ' << ctetraIndex++;
		os << ' ' << squashIndex(get(cell_pmap, cit), squashRanges);
		for (int i = 0; i < 4; i++)
			os << ' ' << V[cit->vertex(i)];
		os << '\n';
	}

	//-------------------------------------------------------
	// Facets
	//-------------------------------------------------------
	typename C3T3::size_type number_of_triangles = c3t3.number_of_facets_in_complex();

	if ( print_each_facet_twice )
		number_of_triangles += number_of_triangles;

	os << "$ Triangles: " << number_of_triangles << '\n';

	int ctriaxIndex = 1;
	for ( Facet_iterator fit = c3t3.facets_in_complex_begin();
			fit != c3t3.facets_in_complex_end();
			++fit)
	{
		typename C3T3::Facet f = (*fit);

		// Apply priority among subdomains, to get consistent facet orientation per subdomain-pair interface.
		if ( print_each_facet_twice )
		{
			// NOTE: We mirror a facet when needed to make it consistent with No_patch_facet_pmap_first/second.
			if (f.first->subdomain_index() > f.first->neighbor(f.second)->subdomain_index())
				f = tr.mirror_facet(f);
		}

		// Get facet vertices in CCW order.
		Vertex_handle vh1 = f.first->vertex((f.second + 1) % 4);
		Vertex_handle vh2 = f.first->vertex((f.second + 2) % 4);
		Vertex_handle vh3 = f.first->vertex((f.second + 3) % 4);

		// Facet orientation also depends on parity.
		if (f.second % 2 != 0)
			std::swap(vh2, vh3);

		os << "CTRIAX";
		os << ' ' << ctriaxIndex++;
		os << ' ' << squashIndex(get(facet_pmap, *fit), squashRanges);
		os << ' ' << V[vh1] << ' ' << V[vh2] << ' ' << V[vh3];
		os << '\n';

		// Print triangle again if needed, with opposite orientation
		if ( print_each_facet_twice )
		{
			os << "CTRIAX";
			os << ' ' << ctriaxIndex++;
			os << ' ' << squashIndex(get(facet_twice_pmap, *fit), squashRanges);
			os << ' ' << V[vh3] << ' ' << V[vh2] << ' ' << V[vh1];
			os << '\n';
		}
	}
}

template <class C3T3,
		class Vertex_index_property_map,
		class Facet_index_property_map,
		class Facet_index_property_map_twice,
		class Cell_index_property_map>
void output_to_nastran(std::ostream & os,
		const C3T3 & c3t3,
		const Vertex_index_property_map & vertex_pmap,
		const Facet_index_property_map & facet_pmap,
		const Cell_index_property_map & cell_pmap,
		const Facet_index_property_map_twice & facet_twice_pmap = Facet_index_property_map_twice(),
		const bool print_each_facet_twice = false)
{
	typedef typename C3T3::Triangulation Tr;
	typedef typename C3T3::Facets_in_complex_iterator Facet_iterator;
	typedef typename C3T3::Cells_in_complex_iterator Cell_iterator;

	typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
	typedef typename Tr::Vertex_handle Vertex_handle;
	typedef typename Tr::Weighted_point Weighted_point;

	const Tr & tr = c3t3.triangulation();

	//-------------------------------------------------------
	// File output
	//-------------------------------------------------------

	//-------------------------------------------------------
	// Header
	//-------------------------------------------------------
	os << std::setprecision(17);

	os << "$ Created by CGAL with MP exporter\n";


	//-------------------------------------------------------
	// Vertices
	//-------------------------------------------------------

	os << "$ Vertices: " << tr.number_of_vertices() << '\n';

	boost::unordered_map<Vertex_handle, int> V;
	int inum = 1;
	for ( Finite_vertices_iterator vit = tr.finite_vertices_begin();
			vit != tr.finite_vertices_end();
			++vit)
	{
		V[vit] = inum++;
		Weighted_point p = tr.point(vit);
		os << "GRID"
				<< ' ' << V[vit]
				<< ' ' << get(vertex_pmap, vit)
				<< ' ' << CGAL::to_double(p.x())
				<< ' ' << CGAL::to_double(p.y())
				<< ' ' << CGAL::to_double(p.z())
				<< '\n';
	}

	//-------------------------------------------------------
	// Tetrahedra
	//-------------------------------------------------------
	os << "$ Tetrahedra: " << c3t3.number_of_cells_in_complex() << '\n';

	int ctetraIndex = 1;
	for ( Cell_iterator cit = c3t3.cells_in_complex_begin() ;
			cit != c3t3.cells_in_complex_end() ;
			++cit )
	{
		os << "CTETRA";
		os << ' ' << ctetraIndex++;
		os << ' ' << get(cell_pmap, cit);
		for (int i = 0; i < 4; i++)
			os << ' ' << V[cit->vertex(i)];
		os << '\n';
	}

	//-------------------------------------------------------
	// Facets
	//-------------------------------------------------------
	typename C3T3::size_type number_of_triangles = c3t3.number_of_facets_in_complex();

	if ( print_each_facet_twice )
		number_of_triangles += number_of_triangles;

	os << "$ Triangles: " << number_of_triangles << '\n';

	int ctriaxIndex = 1;
	for ( Facet_iterator fit = c3t3.facets_in_complex_begin();
			fit != c3t3.facets_in_complex_end();
			++fit)
	{
		typename C3T3::Facet f = (*fit);

		// Apply priority among subdomains, to get consistent facet orientation per subdomain-pair interface.
		if ( print_each_facet_twice )
		{
			// NOTE: We mirror a facet when needed to make it consistent with No_patch_facet_pmap_first/second.
			if (f.first->subdomain_index() > f.first->neighbor(f.second)->subdomain_index())
				f = tr.mirror_facet(f);
		}

		// Get facet vertices in CCW order.
		Vertex_handle vh1 = f.first->vertex((f.second + 1) % 4);
		Vertex_handle vh2 = f.first->vertex((f.second + 2) % 4);
		Vertex_handle vh3 = f.first->vertex((f.second + 3) % 4);

		// Facet orientation also depends on parity.
		if (f.second % 2 != 0)
			std::swap(vh2, vh3);

		os << "CTRIAX";
		os << ' ' << ctriaxIndex++;
		os << ' ' << get(facet_pmap, *fit);
		os << ' ' << V[vh1] << ' ' << V[vh2] << ' ' << V[vh3];
		os << '\n';

		// Print triangle again if needed, with opposite orientation
		if ( print_each_facet_twice )
		{
			os << "CTRIAX";
			os << ' ' << ctriaxIndex++;
			os << ' ' << get(facet_twice_pmap, *fit);
			os << ' ' << V[vh3] << ' ' << V[vh2] << ' ' << V[vh1];
			os << '\n';
		}
	}
}

}

template <class C3T3>
void output_to_nastran(std::ostream & os,
		const C3T3 & c3t3,
		bool rebind = false,
		bool show_patches = false)
{
	if (rebind) {
		if (show_patches)
			Mesh_3::output_to_nastran<C3T3, true, false>(os, c3t3);
		else
			Mesh_3::output_to_nastran<C3T3, true, true>(os, c3t3);
	} else {
		if (show_patches)
			Mesh_3::output_to_nastran<C3T3, false, false>(os, c3t3);
		else
			Mesh_3::output_to_nastran<C3T3, false, true>(os, c3t3);
	}
}

template <class C3T3>
void squashed_output_to_nastran(std::ostream & os,
		const C3T3 & c3t3,
		bool rebind = false,
		bool show_patches = false,
		const std::vector<std::pair<int, int>> & squashRanges = {})
{
	if (rebind) {
		if (show_patches)
			Mesh_3::squashed_output_to_nastran<C3T3, true, false>(os, c3t3, squashRanges);
		else
			Mesh_3::squashed_output_to_nastran<C3T3, true, true>(os, c3t3, squashRanges);
	} else {
		if (show_patches)
			Mesh_3::squashed_output_to_nastran<C3T3, false, false>(os, c3t3, squashRanges);
		else
			Mesh_3::squashed_output_to_nastran<C3T3, false, true>(os, c3t3, squashRanges);
	}
}

}

#endif // NASTRAN_EXPORT_HPP
