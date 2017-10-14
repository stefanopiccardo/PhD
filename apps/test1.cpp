#include <iostream>
#include "src/mesh.hpp"
#include "src/silo_io.hpp"

int main(int argc, char **argv)
{
    using   coord_type = double;
    using   scalar_type = double;

    if (argc != 2)
    {
        std::cout << "Usage: " << argv[0] << " <mesh file name>" << std::endl;
        return 1;
    }

    mesh<coord_type> msh;
    load_netgen_2d(argv[1], msh);

    silo_database silo;
    silo.create("test.silo");
    silo.add_mesh(msh, "test_mesh");


    auto level_set_function = [](const typename mesh<coord_type>::point_type pt) -> coord_type {
        auto x = pt.x();
        auto y = pt.y();
        auto alpha = 0.4;
        auto beta = 0.4;

        return (x-alpha)*(x-alpha) + (y-beta)*(y-beta) - 0.1;
    };

    std::vector<coord_type> vals;

    for (auto& node : msh.nodes)
    {
        auto pt = points(msh, node);
        auto val = level_set_function(pt);
        vals.push_back(val);
    }

    silo.add_variable("test_mesh", "level_set", vals.data(), vals.size(), nodal_variable_t);

    detect_cut_cells(msh, level_set_function);

    std::vector<coord_type> marked_cells;
    marked_cells.reserve( msh.cells.size() );

    for (auto& cl : msh.cells)
        if (cl.cut_by_interface)
            marked_cells.push_back(1.0);
        else
            marked_cells.push_back(0.0);

    silo.add_variable("test_mesh", "cut", marked_cells.data(), marked_cells.size(), zonal_variable_t);

    silo.close();

    return 0;
}
