#ifndef UTILS_H
#define UTILS_H
#include <zlib.h>

typedef struct String String;
struct String
{
	size_t size;
	size_t capacity;
	char* str;
};

String* string_create(void);

void string_destroy(String* str);
void string_reset(String* str);

char* parse_line(gzFile fh, String* str);

void line_destroy();
#endif

