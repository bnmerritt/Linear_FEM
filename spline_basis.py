import numpy
import matplotlib
import readBEXT_JSON
import basis

def evaluateElementBernsteinBasisAtParamCoord( uspline, elem_id, param_coord ):
    elem_degree = readBEXT_JSON.getElementDegree(uspline,elem_id) # Get the degree of the element
    elem_bernstein_basis = numpy.zeros( elem_degree + 1 )
    elem_domain = readBEXT_JSON.getElementDomain( uspline, elem_id )
    for n in range( 0, elem_degree + 1 ):
        elem_bernstein_basis[n] = basis.evalBernsteinBasis1D(param_coord, elem_degree, n, elem_domain)# Evaluate the Bernstein basis at the parametric coordinate
    # print(elem_bernstein_basis)
    return elem_bernstein_basis

def evaluateElementSplineBasisAtParamCoord( uspline, elem_id, param_coord ):
    elem_ext_operator = readBEXT_JSON.getElementExtractionOperator( uspline, elem_id ) # Get the extraction operator of the element
    elem_bernstein_basis = evaluateElementBernsteinBasisAtParamCoord( uspline, elem_id, param_coord )
    # Main issue is how to solve this
    elem_spline_basis = sum(numpy.transpose(elem_ext_operator * elem_bernstein_basis)) # Apply the extraction operator onto its Bernstein basis at the param coord
    return elem_spline_basis 

def plotUsplineBasis( uspline, color_by ):
    num_pts = 100
    xi = numpy.linspace( 0, 1, num_pts )
    num_elems = readBEXT_JSON.getNumElems( uspline ) # Get number of elements in the U-spline
    for elem_idx in range( 0, num_elems ):
        elem_id = readBEXT_JSON.elemIdFromElemIdx( uspline, elem_idx ) # Get the element id of the current element
        elem_domain = readBEXT_JSON.getElementDomain( uspline, elem_id ) # Get the domain of the current element
        elem_node_ids = readBEXT_JSON.getElementNodeIds( uspline, elem_id ) # Get the spline node ids of the current element
        elem_degree = readBEXT_JSON.getElementDegree(uspline,elem_id) # I added this because y needed the degree
        x = numpy.linspace( elem_domain[0], elem_domain[1], num_pts )
        y = numpy.zeros( shape = ( elem_degree + 1, num_pts ) )
        for i in range( 0, num_pts ):
            y[:,i] = evaluateElementSplineBasisAtParamCoord( uspline, elem_id, x) # Evaluate the current element's spline basis at the current coordinate
        # Do plotting for the current element
        for n in range( 0, elem_degree + 1 ):
            if color_by == "element":
                color = getLineColor( elem_idx )
            elif color_by == "node":
                color = getLineColor( elem_node_ids[n] )
            matplotlib.pyplot.plot( x, y[n,:], color = getLineColor( color ) )
    matplotlib.plt.show()

def getLineColor( idx ):
    colors = list( matplotlib.colors.TABLEAU_COLORS.keys() )
    num_colors = len( colors )
    mod_idx = idx % num_colors
    return matplotlib.colors.TABLEAU_COLORS[ colors[ mod_idx ] ]

uspline = readBEXT_JSON.readBEXT( "quadratic_bspline_match.json" )
num_pts = 100
color_by = "element"
xi = numpy.linspace( 0, 1, num_pts )
num_elems = readBEXT_JSON.getNumElems( uspline ) # Get number of elements in the U-spline
for elem_idx in range( 0, num_elems):
    elem_id = readBEXT_JSON.elemIdFromElemIdx( uspline, elem_idx ) # Get the element id of the current element
    elem_domain = readBEXT_JSON.getElementDomain( uspline, elem_id ) # Get the domain of the current element
    elem_node_ids = readBEXT_JSON.getElementNodeIds( uspline, elem_id ) # Get the spline node ids of the current element
    elem_degree = readBEXT_JSON.getElementDegree(uspline,elem_id) # I added this because y needed the degree
    x = numpy.linspace( elem_domain[0], elem_domain[1], num_pts )
    print(elem_domain)
    y = numpy.zeros( shape = ( elem_degree + 1, num_pts ) )
    for i in range( 0, num_pts ):
        y[:,i] = evaluateElementSplineBasisAtParamCoord( uspline, elem_id, x[i]) # Evaluate the current element's spline basis at the current coordinate
    # Do plotting for the current element
    for n in range( 0, elem_degree + 1 ):
        if color_by == "element":
            color = getLineColor( elem_idx )
        elif color_by == "node":
            color = getLineColor( elem_node_ids[n] )
        matplotlib.pyplot.plot( x, y[n,:])
# matplotlib.pyplot.plt.show()

# uspline = readBEXT_JSON.readBEXT( "quadratic_bspline_match.json" )
# plotUsplineBasis( uspline, "element" )
# plotUsplineBasis( uspline, "node" )

# Check to make sure I am summing in rows not columns
# A = numpy.asarray([[1,0,0],
#       [0,1,0.5],
#       [0,0,0.5]])
# B = numpy.asarray([[0.25],[0.5],[0.25]])
# C1 = A * B
# C = sum(numpy.transpose(A * B))

# uspline = readBEXT_JSON.readBEXT( "quadratic_bspline_match.json" )
# elem_id = 0
# num_pts = 100
# param_coord = numpy.linspace(-1,1,num_pts)

# for elem_idx in range( 0, 1):
#     elem_id = readBEXT_JSON.elemIdFromElemIdx( uspline, elem_idx ) # Get the element id of the current element
#     elem_domain = readBEXT_JSON.getElementDomain( uspline, elem_id ) # Get the domain of the current element
#     elem_node_ids = readBEXT_JSON.getElementNodeIds( uspline, elem_id ) # Get the spline node ids of the current element
#     elem_degree = readBEXT_JSON.getElementDegree(uspline,elem_id) # I added this because y needed the degree
#     x = numpy.linspace( elem_domain[0], elem_domain[1], num_pts )
#     print(elem_domain)
#     y = numpy.zeros( shape = ( elem_degree + 1, num_pts ) )
#     for i in range( 0, num_pts ):
#         y[:,i] = evaluateElementSplineBasisAtParamCoord( uspline, elem_id, x[i]) # Evaluate the current element's spline basis at the current coordinate

# elem_ext_operator = readBEXT_JSON.getElementExtractionOperator( uspline, elem_id ) # Get the extraction operator of the element
# elem_bernstein_basis = evaluateElementBernsteinBasisAtParamCoord( uspline, elem_id, param_coord )
# elem_spline_basis = numpy.multiply(elem_ext_operator,elem_bernstein_basis) # Apply the extraction operator onto its Bernstein basis at the param coord