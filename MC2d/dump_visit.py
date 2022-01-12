import numpy as np

def dump_visit(infile, points, triangles):
    """
    Snippet to write data of a mesh in a specific
    format to be able to load in visit 
    Download visit binaries from
    https://wci.llnl.gov/simulation/computer-codes/visit/executables
    To visualize
    load <output>.vtk in visit
    in Add select subsets>domains
    click draw in gui
    """
    Np = np.shape(points)[0]
    num_triangles = np.shape(triangles)[0]
    print (Np, num_triangles) 
    ## write the headers of the file
    with open(infile, "w") as f:
        f.write('# vtk DataFile Version 2.0 \n')
        f.write('grid, time 110\n')
        f.write('ASCII \n')
        f.write('DATASET POLYDATA \n')
        f.write('POINTS  '+str(Np)+'  float\n')
        for p in points:
            f.write("%16.8f %16.8f %16.8f\n" %(p[0], p[1], p[2]))
        f.write('POLYGONS  '+str(num_triangles)+'  '
                +str(4*num_triangles) + '\n')
        for it, tri in enumerate(triangles):
            f.write("%d %d %d %d\n" %(3, tri[0], tri[1], tri[2]))

    
