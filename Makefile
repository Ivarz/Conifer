conifer : src/main.c src/utils.c src/kraken_stats.c src/kraken_taxo.c
	gcc -std=c99 -Wall -Wextra -O3 -D_POSIX_C_SOURCE=200809L -I third_party/uthash/src -I . src/utils.c src/kraken_stats.c src/kraken_taxo.c src/main.c -o conifer -lm -lz

tests : test/tests.c src/kraken_stats.c src/kraken_taxo.c
	cmake third_party/Unity/CMakeLists.txt
	make -C third_party/Unity/
	gcc -std=c99 -Wall -D_POSIX_C_SOURCE=200809L -I third_party/uthash/src -I . -I third_party/Unity/src src/kraken_stats.c src/kraken_taxo.c src/utils.c $< -L third_party/Unity/src -o tests -lunity -lm -lz
	./tests; rm tests

is_a_parent_of_b : utils/is_a_parent_of_b.c src/kraken_stats.c src/kraken_taxo.c
	gcc -std=c99 -Wall -O3 -D_POSIX_C_SOURCE=200809L -I third_party/uthash/src -I . src/kraken_stats.c src/kraken_taxo.c utils/is_a_parent_of_b.c -o is_a_parent_of_b -lm

taxid_name : utils/taxid_name.c src/kraken_stats.c src/kraken_taxo.c
	gcc -std=c99 -Wall -O3 -D_POSIX_C_SOURCE=200809L -I third_party/uthash/src -I . src/kraken_stats.c src/kraken_taxo.c utils/taxid_name.c -o taxid_name -lm

show_ancestors : utils/show_ancestors.c src/kraken_stats.c src/kraken_taxo.c
	gcc -std=c99 -Wall -O3 -D_POSIX_C_SOURCE=200809L -I third_party/uthash/src -I . src/kraken_stats.c src/kraken_taxo.c utils/show_ancestors.c -o show_ancestors -lm
