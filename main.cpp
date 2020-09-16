#define CGAL_CONCURRENT_MESH_3 1

#include "crap.hpp"

#include "nastran_export.hpp"

#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/helpers.h>
//#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>

#include <cstdlib>
#include <filesystem>
#include <string>

struct PolyDomainArgs
{
	// Input output.
	std::string inputFile;
	std::string meditFile;
	std::string bdfFile;

	// Mesh criteria.
	double edgeSize;
	double facetAngle;
	double facetSize;
	double facetDistance;
	double cellRadiusEdgeRatio;
	double cellSize;
};

#ifdef CGAL_CONCURRENT_MESH_3
typedef CGAL::Parallel_tag ConcurrencyTag;
#else
typedef CGAL::Sequential_tag ConcurrencyTag;
#endif

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Mesh_polyhedron_3<K>::type Polyhedron;

void readPolyhedron(Polyhedron & polyhedron, const std::string & filePath)
{
	std::ifstream input(filePath);
	input >> polyhedron;
	if (input.fail()) {
		std::cerr << "Error: Could not read file '" <<  filePath << "'" << std::endl;
		std::exit(EXIT_FAILURE);
	}
	input.close();

	CGAL::Polygon_mesh_processing::triangulate_faces(polyhedron);

	if (!CGAL::is_triangle_mesh(polyhedron)) {
		std::cerr << "Input geometry is not triangulated." << std::endl;
		std::exit(EXIT_FAILURE);
	}
}

void polyDomain(const PolyDomainArgs & args)
{
	// Domain
	typedef CGAL::Polyhedral_mesh_domain_with_features_3<K> Domain;

	// Triangulation
	typedef CGAL::Mesh_triangulation_3<Domain, CGAL::Default, ConcurrencyTag>::type Tr;
	typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3T3;

	// Criteria
	typedef CGAL::Mesh_criteria_3<Tr> MeshCriteria;

	// To avoid verbose function and named parameters call
	using namespace CGAL::parameters;

	Polyhedron casePoly;
	readPolyhedron(casePoly, args.inputFile);

	// Create domain
	Domain domain(casePoly);
	domain.detect_features();

	MeshCriteria criteria(edge_size = args.edgeSize,
			facet_angle = args.facetAngle,
			facet_size = args.facetSize,
			facet_distance = args.facetDistance,
			cell_radius_edge_ratio = args.cellRadiusEdgeRatio,
			cell_size = args.cellSize);

	// Mesh generation
	C3T3 c3t3 = CGAL::make_mesh_3<C3T3>(domain, criteria, no_perturb(), no_exude());

	// Output
	std::ofstream medit_file(args.meditFile);
	c3t3.output_to_medit(medit_file, false, true);
	medit_file.close();

	// Output nastran
	std::ofstream nastran_file(args.bdfFile);
	CGAL::output_to_nastran(nastran_file, c3t3, false, true);
	medit_file.close();
}

int main(int argc, char * argv[])
{
	crap::KeyArg programArg(argv[0]);
	programArg.setRequired(true);

	crap::Parser parser(& programArg);
	parser.setOptionRequired(true);
	parser.setHeader("CGAL based BDF mesh generator\n\n");
	parser.setFooter("\n(c)2020 Michal Policht. WTFPL.\n");

	crap::KeyArg helpArg("help", "Print this information");
	helpArg.addAlias("--help").addAlias("-h");
	parser.addSubCmd(& helpArg);

	crap::KeyArg generateArg("generate", "Generate mesh");
	crap::ArgGroup generateArgsGroup("options");
	crap::ValueArg inputFileArg("input_file", "Input file.");
	inputFileArg.setRequired(true);
	generateArgsGroup.addAttr(& inputFileArg);
	crap::ValueArg outputDirArg("output_dir", "Output directory.");
	outputDirArg.setDefaultValue(".");
	generateArgsGroup.addAttr(& outputDirArg);
	crap::KeyValueArg edgeSizeArg("edge_size", "size", "Edge <size>.");
	edgeSizeArg.setDefaultValue("0.5");
	generateArgsGroup.addAttr(& edgeSizeArg);
	crap::KeyValueArg facetSizeArg("facet_size", "size", "Facet <size>.");
	facetSizeArg.setDefaultValue("0.5");
	generateArgsGroup.addAttr(& facetSizeArg);
	crap::KeyValueArg facetAngleArg("facet_angle", "angle", "Facet <angle>.");
	facetAngleArg.setDefaultValue("25");
	generateArgsGroup.addAttr(& facetAngleArg);
	crap::KeyValueArg facetDistanceArg("facet_distance", "distance", "Facet <distance>.");
	facetDistanceArg.setDefaultValue("0.125");
	generateArgsGroup.addAttr(& facetDistanceArg);
	crap::KeyValueArg cellRadiusEdgeRatioArg("cell_radius_edge_ratio", "ratio", "Cell radius edge <ratio>.");
	cellRadiusEdgeRatioArg.setDefaultValue("3");
	generateArgsGroup.addAttr(& cellRadiusEdgeRatioArg);
	crap::KeyValueArg cellSizeArg("cell_size", "size", "Cell <size>.");
	cellSizeArg.setDefaultValue("2.5");
	generateArgsGroup.addAttr(& cellSizeArg);
	parser.addSubCmd(& generateArg)->addArgGroup(& generateArgsGroup);

	try {
		parser.parse(argc, argv);
	} catch (const crap::Exception & e) {
		std::cout << "\n" << e.what() << "\n\n";
		parser.printSynopsis();
		std::cout << "\n";
		return EXIT_FAILURE;
	}

	if (helpArg.isSet())
		parser.printHelp(std::cout);
	else if (generateArg.isSet()) {
		PolyDomainArgs args;
		args.inputFile = inputFileArg.value();

		std::string outputDir(".");
		if (outputDirArg.isSet())
			outputDir = outputDirArg.value();
		args.meditFile = (std::filesystem::path(outputDir) / (std::filesystem::path(args.inputFile).stem().string() + ".mesh")).string();
		args.bdfFile = (std::filesystem::path(outputDir) / (std::filesystem::path(args.inputFile).stem().string() + ".bdf")).string();

		args.edgeSize = std::stod(edgeSizeArg.value());
		args.facetAngle = std::stod(facetAngleArg.value());
		args.facetSize = std::stod(facetSizeArg.value());
		args.facetDistance = std::stod(facetDistanceArg.value());
		args.cellRadiusEdgeRatio = std::stod(cellRadiusEdgeRatioArg.value());
		args.cellSize = std::stod(cellSizeArg.value());

		std::cout << "Input file: " << args.inputFile << "\n";
		std::cout << "Output directory: " << outputDir << "\n";
		std::cout << "Medit output file: " << args.meditFile << "\n";
		std::cout << "Nastran output file: " << args.bdfFile << "\n";
		std::cout << "Edge size: " << args.edgeSize << "\n";
		std::cout << "Facet angle: " << args.facetAngle << "\n";
		std::cout << "Facet size: " << args.facetSize << "\n";
		std::cout << "Facet distance: " << args.facetDistance << "\n";
		std::cout << "Cell radius edge ratio: " << args.cellRadiusEdgeRatio << "\n";
		std::cout << "Cell size: " << args.cellSize << "\n";

		polyDomain(args);
	}

	return EXIT_SUCCESS;
}
