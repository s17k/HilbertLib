struct MDPoint {
	coord_t* coordinates;
	id_t own_data_id;
};


typedef struct MDPoint MDPoint;

void make_MDPoint(MDPoint *X, int dimensions) {
	X->coordinates = calloc(dimensions, sizeof(coord_t));
	X->own_data_id = 0;
};

void MDPointRemove(MDPoint *X) {
	free(X->coordinates);
	X->own_data_id = 0;
}


