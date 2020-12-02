#include <stdlib.h>
#include "src/utils.h"
#include <string.h>
#include <stdio.h>

String* string_create(void)
{     
	String* str = malloc(sizeof(*str));
    char* content = malloc(sizeof(*content)*4096);
	str->capacity = 4096;
	str->size = 0;
	str->str = content;
	memset(str->str, 0, sizeof(*str->str)*str->capacity);
    return str;
}

String* string_create_from(char* s)
{
	size_t size = strnlen(s, 4096);

	String* str = string_create();
	str->size = size;
	str->str = strncpy(str->str, s, sizeof(*str->str)*size);
	if (str->str[size] != '\0'){
		fprintf(stderr, "WARNING! Creating String from char* longer than buffersize of 4095\n");
	}
	return str;
}

String* string_create_copy(String* str)
{     
	String* cpy = malloc(sizeof(*cpy));
	cpy->capacity = str->capacity;
	cpy->size = str->size;
	cpy->str = malloc(sizeof(*str->str)*str->capacity);
	cpy->str = strncpy(cpy->str, str->str, sizeof(*str->str)*str->capacity);
    return cpy;
}

void string_destroy(String* s)
{
	s->size = 0;	
	free(s->str);
	free(s);
}

void string_reset(String* s)
{
	memset(s->str, 0, sizeof(*s->str)*(s->size+1));
	return;
}

void remove_trailing_newline(String* s)
{
	size_t last_idx = s->size-1;
	if (s->str[last_idx] == '\n'){
		s->str[last_idx] = '\0';
		s->size--;
	}
	return;
}
char* parse_line(gzFile fh, String* s)
{
	char* status = gzgets(fh, s->str, sizeof(*s->str)*s->capacity);
	s->size = strnlen(s->str, s->capacity);
	remove_trailing_newline(s);
	while (s->size + 1 == s->capacity){
		s->capacity = 2*(s->capacity);
		s->str = realloc(s->str, sizeof(*s->str)*s->capacity);
		status = gzgets(fh, s->str + s->size, sizeof(*s->str)*(s->capacity - s->size));
		s->size = strnlen(s->str, s->capacity);
		remove_trailing_newline(s);
	}
	return status;
}
