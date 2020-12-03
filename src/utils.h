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
String* string_create_from(char* str);
String* string_create_copy(String* str);
void string_copy(String* dst, String const* const src);

void string_destroy(String* str);
void string_reset(String* str);

char* parse_line(gzFile fh, String* str);

void line_destroy();
#endif

