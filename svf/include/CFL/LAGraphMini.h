typedef enum    // GrB_Info
{

    GrB_SUCCESS = 0,            // all is well

    //--------------------------------------------------------------------------
    // informational codes, not an error:
    //--------------------------------------------------------------------------

    GrB_NO_VALUE = 1,           // A(i,j) requested but not there
    #ifndef GRAPHBLAS_VANILLA
    GxB_EXHAUSTED = 7089,       // iterator is exhausted
    #endif

    //--------------------------------------------------------------------------
    // errors:
    //--------------------------------------------------------------------------

    GrB_UNINITIALIZED_OBJECT = -1,  // object has not been initialized
    GrB_NULL_POINTER = -2,          // input pointer is NULL
    GrB_INVALID_VALUE = -3,         // general error; some value is bad
    GrB_INVALID_INDEX = -4,         // row or column index is out of bounds
    GrB_DOMAIN_MISMATCH = -5,       // object domains are not compatible
    GrB_DIMENSION_MISMATCH = -6,    // matrix dimensions do not match
    GrB_OUTPUT_NOT_EMPTY = -7,      // output matrix already has values
    GrB_NOT_IMPLEMENTED = -8,       // method not implemented
    GrB_ALREADY_SET = -9,           // field already written to
    GrB_PANIC = -101,               // unknown error
    GrB_OUT_OF_MEMORY = -102,       // out of memory
    GrB_INSUFFICIENT_SPACE = -103,  // output array not large enough
    GrB_INVALID_OBJECT = -104,      // object is corrupted
    GrB_INDEX_OUT_OF_BOUNDS = -105, // row or col index out of bounds
    GrB_EMPTY_OBJECT = -106,        // an object does not contain a value
    #ifndef GRAPHBLAS_VANILLA
    GxB_JIT_ERROR = -7001,          // JIT compiler/loader error
    GxB_GPU_ERROR = -7002,          // GPU error (future; not yet in production)
    GxB_OUTPUT_IS_READONLY = -7003, // output matrix has readonly components
    #endif

}
GrB_Info ;

typedef uint64_t GrB_Index ;
typedef struct GB_Matrix_opaque       *GrB_Matrix ;
typedef struct GB_Type_opaque         *GrB_Type ;
typedef struct GB_BinaryOp_opaque     *GrB_BinaryOp ;

 typedef struct {
    int32_t nonterm; // prod_A != -1 && prod_B != -1 => Type of Rule is [Variable -> AB]
    int32_t prod_A;  // prod_A == -1 && prod_B == -1 => Type of Rule is [Variable -> eps]
    int32_t prod_B;  // prod_A != -1 && prod_B == -1 => Type of Rule is [Variable -> term]
    int32_t index;   // For rules that can be grouped by index
 } LAGraph_rule_WCNF;

    GrB_Info GrB_Matrix_new     // create a new matrix with no entries
(
    GrB_Matrix *A,          // handle of matrix to create
    GrB_Type type,          // type of matrix to create
    GrB_Index nrows,        // matrix dimension is nrows-by-ncols
    GrB_Index ncols         // (nrows and ncols must be <= GrB_INDEX_MAX+1)
) ;

 GrB_Info GrB_Matrix_build  // build a matrix from (I,J,X) tuples
  (
      GrB_Matrix C,               // matrix to build
      const GrB_Index *I,         // array of row indices of tuples
      const GrB_Index *J,         // array of column indices of tuples
      const bool *X,            // array of values of tuples
      GrB_Index nvals,            // number of tuples
      const GrB_BinaryOp dup      // binary function to assemble duplicates
  ) ;

GrB_Info GrB_Matrix_nvals   // get the number of entries in a matrix
(
    GrB_Index *nvals,       // matrix has nvals entries
    const GrB_Matrix A      // matrix to query
) ;

GrB_Info GrB_Matrix_extractTuples     // [I,J,X] = find (A)
(
      uint64_t *I,            // array for returning row indices of tuples
      uint64_t *J,            // array for returning col indices of tuples
      bool *X,              // array for returning values of tuples
      GrB_Index *nvals,       // I,J,X size on input; # tuples on output
      const GrB_Matrix A      // matrix to extract tuples from
  ) ;

GrB_Info LAGraph_CFL_reachability
(
    // Output
    GrB_Matrix *outputs, // Array of matrices containing results.
                         // The size of the array must be equal to nonterms_count.
                         //
                         // outputs[k]: (i, j) = true if and only if there is a path
                         // from node i to node j whose edge labels form a word
                         // derivable from the non-terminal 'k' of the specified CFG.
    // Input
    const GrB_Matrix *adj_matrices, // Array of adjacency matrices representing the graph.
                                    // The length of this array is equal to the count of
                                    // terminals (terms_count).
                                    //
                                    // adj_matrices[t]: (i, j) == 1 if and only if there
                                    // is an edge between nodes i and j with the label of
                                    // the terminal corresponding to index 't' (where t is
                                    // in the range [0, terms_count - 1]).
    int64_t terms_count,            // The total number of terminal symbols in the CFG.
    int64_t nonterms_count,         // The total number of non-terminal symbols in the CFG.
    const LAGraph_rule_WCNF *rules, // The rules of the CFG.
    int64_t rules_count,            // The total number of rules in the CFG.
    char *msg                       // Message string for error reporting.
) ;
