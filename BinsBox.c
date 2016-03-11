struct BinsBox
{
	coord_t *data;
	size_t size;
	size_t capacity;
};

typedef struct BinsBox BinsBox;

void
BinsBoxDoubleCapacity (BinsBox * b)
{
	b->data =
		(coord_t *) realloc ((void *) b->data,
				     b->capacity * 2);
	b->capacity *= 2;
}

void
BinsBoxPushBack (BinsBox * b, coord_t x)
{
	b->size++;
	if (b->size > b->capacity)
		BinsBoxDoubleCapacity (b);
	b->data[b->size - 1] = x;
}

coord_t
BinsBoxGet (BinsBox * b, size_t index)
{
	return b->data[index];
}

void
BinsBoxDeallocate (BinsBox * b)
{
	free (b->data);
	b->size = 0;
	b->capacity = 0;
}

typedef struct BinsBox BinsBox;
