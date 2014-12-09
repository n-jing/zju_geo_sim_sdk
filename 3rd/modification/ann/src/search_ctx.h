#ifndef HJ_ANN_SEARCH_CTX_H_
#define HJ_ANN_SEARCH_CTX_H_

class ANNmin_k;
class ANNPr_queue;

class search_ctx {
public:
	search_ctx()
		:ANNmaxPtsVisited (0){}
	int	ANNmaxPtsVisited;	// maximum number of pts visited
	int	ANNptsVisited;			// number of pts visited in search
};
class annkFRSearch_ctx : public search_ctx {
public:
	int				ANNkdFRDim;				// dimension of space
	ANNpoint		ANNkdFRQ;				// query point
	ANNdist			ANNkdFRSqRad;			// squared radius search bound
	double			ANNkdFRMaxErr;			// max tolerable squared error
	ANNpointArray	ANNkdFRPts;				// the points
	ANNmin_k*		ANNkdFRPointMK;			// set of k closest points
	int				ANNkdFRPtsVisited;		// total points visited
	int				ANNkdFRPtsInRange;		// number of points in the range
};
class ann_search_ctx : public search_ctx {
public:
	int				ANNkdDim;				// dimension of space
	ANNpoint		ANNkdQ;					// query point
	double			ANNkdMaxErr;			// max tolerable squared error
	ANNpointArray	ANNkdPts;				// the points
	ANNmin_k		*ANNkdPointMK;			// set of k closest points
};
class annkPriSearch_ctx : public search_ctx {
public:
	double			ANNPrEps;				// the error bound
	int				ANNPrDim;				// dimension of space
	ANNpoint		ANNPrQ;					// query point
	double			ANNPrMaxErr;			// max tolerable squared error
	ANNpointArray	ANNPrPts;				// the points
	ANNPr_queue		*ANNPrBoxPQ;			// priority queue for boxes
	ANNmin_k		*ANNPrPointMK;			// set of k closest points
};

#endif
