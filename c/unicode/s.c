#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <windows.h>
#define _root_pai 3.141592653589793
#define _root_e 2.718281828459045
#define _root_small 1.0e-12
#define _root_large 1.0e12
typedef void _root_void;
typedef char _root_character_1;
typedef int8_t _root_integer_1;
typedef int16_t _root_integer_2;
typedef int32_t _root_integer_4;
typedef int64_t _root_integer_8;
typedef uint8_t _root_integer_1_unsigned;
typedef uint16_t _root_integer_2_unsigned;
typedef uint32_t _root_integer_4_unsigned;
typedef uint64_t _root_integer_8_unsigned;
typedef float _root_float_4;
typedef double _root_float_8;
_root_void _root_convert_utf16_utf8(_root_integer_2_unsigned *_to, _root_integer_4 _to_count, _root_integer_1_unsigned *_from) {
	_root_integer_4 _to_index = 0, _from_index = 0;
	while (_from[_from_index] && _to_index < _to_count) {
		_root_integer_4 _code = _from[_from_index++];
		if (_code < 0x80) {
			_to[_to_index++] = _code;
			continue;
		}
		if (_code < 0xc0) {
			++_from_index;
			continue;
		}
		if (_code < 0xe0) {
			_root_integer_4 _to_code = (_code & 0x1f) << 6;
			_to_code |= _from[_from_index++] & 0x3f;
			_to[_to_index++] = _to_code;
			continue;
		}
		if (_code < 0xf0) {
			_root_integer_4 _to_code = (_code & 0x0f) << 12;
			_to_code |= (_from[_from_index++] & 0x3f) << 6;
			_to_code |= _from[_from_index++] & 0x3f;
			_to[_to_index++] = _to_code;
			continue;
		}
		if (_code < 0xf8) {
			_root_integer_4 _to_code = (_code & 0x07) << 18;
			_to_code |= (_from[_from_index++] & 0x3f) << 12;
			_to_code |= (_from[_from_index++] & 0x3f) << 6;
			_to_code |= _from[_from_index++] & 0x3f;
			_to[_to_index++] = 0xd800 | _to_code - 0x10000 >> 10;
			_to[_to_index++] = 0xdc00 | _to_code - 0x10000 & 0x3ff;
			continue;
		}
		++_from_index;
	}
	_to[_to_index] = 0;
	return;
}
_root_integer_4 _root_main(_root_integer_4 _arguments_count, _root_character_1 **_arguments) {
	_root_integer_1_unsigned str[] = "Hello World";
	for (int i = 0; str[i]; ++i) {
		printf("%x ", str[i]);
	}
	puts("\n");
	_root_integer_1_unsigned str1[] = u8"あ　亜";
	for (int i = 0; str1[i]; ++i) {
		printf("%x ", str1[i]);
	}
	puts("\n");
	_root_integer_2_unsigned str2[100];
	_root_convert_utf16_utf8(str2, 100, str1);
	for (int i = 0; str2[i]; ++i) {
		printf("%x ", str2[i]);
	}
	puts("\n");
	printf("%ls", L"あ　亜");
	return 0;
}
_root_integer_4 main(_root_integer_4 _arguments_count, _root_character_1 **_arguments) {
	return _root_main(_arguments_count, _arguments);
}
