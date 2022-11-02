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
typedef char _root_character;
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
typedef struct _root_integer_4_2 _root_integer_4_2;
struct _root_integer_4_2 {
	_root_integer_4 _0, _1;
};
_root_void _root_click(_root_integer_4 _right) {
	INPUT _inputs[] = {{INPUT_MOUSE, .mi = {0, 0, 0, _right ? MOUSEEVENTF_RIGHTDOWN : MOUSEEVENTF_LEFTDOWN, 0, 0}}, {INPUT_MOUSE, .mi = {0, 0, 0, _right ? MOUSEEVENTF_RIGHTUP : MOUSEEVENTF_LEFTUP, 0, 0}}};
	SendInput(2, _inputs, sizeof(INPUT));
	return;
}
_root_void _root_type(_root_integer_4 _virtual_key) {
	INPUT _inputs[] = {{INPUT_KEYBOARD, .ki = {_virtual_key, 0, 0, 0, 0}}, {INPUT_KEYBOARD, .ki = {_virtual_key, 0, KEYEVENTF_KEYUP, 0, 0}}};
	SendInput(2, _inputs, sizeof(INPUT));
	return;
}
_root_integer_4 _root_key_state(_root_integer_4 _virtual_key) {
	return GetAsyncKeyState(_virtual_key) & 0x0001;
}
_root_void _root_cursor_position_set(_root_integer_4_2 _position) {
	SetCursorPos(_position._0, _position._1);
	return;
}
_root_integer_4_2 _root_cursor_position_get(void) {
	POINT _position;
	GetCursorPos(&_position);
	return (_root_integer_4_2){_position.x, _position.y};
}
_root_integer_4_2 _root_integer_4_2_add(_root_integer_4_2 _0, _root_integer_4_2 _1) {
	return (_root_integer_4_2){_0._0 + _1._0, _0._1 + _1._1};
}
_root_void _root_window_activate(_root_integer_2_unsigned *_name) {
	HWND _handle = FindWindowW(NULL, _name);
	if (IsIconic(_handle)) {
		ShowWindow(_handle, 9);
		return;
	}
	SetWindowPos(_handle, 0, 0, 0, 0, 0, SWP_NOMOVE | SWP_NOSIZE | SWP_SHOWWINDOW);
	return;
}
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
_root_integer_4 _root_string_count(_root_integer_1_unsigned *_string) {
	return strlen((_root_character *)_string);
}
_root_integer_4 _root_string_find_last(_root_integer_1_unsigned *_string, _root_integer_1_unsigned _key) {
	for (_root_integer_4 _index = _root_string_count(_string) - 1; _index >= 0; --_index) {
		if (_string[_index] == _key) {
			return _index;
		}
	}
	return -1;
}
_root_void _root_string_replace(_root_integer_4 _string_count, _root_integer_1_unsigned *_string, _root_integer_1_unsigned _key, _root_integer_1_unsigned _value) {
	for (_root_integer_4 _index = 0; _index < _string_count; ++_index) {
		if (_string[_index] == _key) {
			_string[_index] = _value;
		}
	}
	return;
}
_root_integer_4 _root_main(_root_integer_4 _parameters_count, _root_character **_parameters) {
	_root_integer_4 _mode = 0;
	_root_integer_4 _row_height = 14;
	_root_integer_4 _row_count = 57;
	_root_integer_4 _row_index = _row_count / 2;
	_root_integer_4 _delay = 1;
	_root_integer_4 _virtual_key_default = 0x88;
	_root_integer_4 _virtual_key_mode_0 = _virtual_key_default;
	_root_integer_4 _virtual_key_mode_1 = _virtual_key_default;
	_root_integer_4 _virtual_key_mode_2 = _virtual_key_default;
	_root_integer_4 _virtual_key_mode_3 = _virtual_key_default;
	_root_integer_4 _virtual_key_price_up = _virtual_key_default;
	_root_integer_4 _virtual_key_price_down = _virtual_key_default;
	_root_integer_4 _virtual_key_buy = _virtual_key_default;
	_root_integer_4 _virtual_key_buy_cancel = _virtual_key_default;
	_root_integer_4 _virtual_key_sell = _virtual_key_default;
	_root_integer_4 _virtual_key_sell_cancel = _virtual_key_default;
	_root_integer_4 _virtual_key_scroll_up = _virtual_key_default;
	_root_integer_4 _virtual_key_scroll_down = _virtual_key_default;
	_root_integer_4 _virtual_key_entry = _virtual_key_default;
	_root_integer_4 _virtual_key_exit = _virtual_key_default;
	_root_integer_4 _virtual_key_volume_up = _virtual_key_default;
	_root_integer_4 _virtual_key_volume_down = _virtual_key_default;
	_root_integer_4 _virtual_key_enter = _virtual_key_default;
	_root_integer_4 _virtual_key_escape = _virtual_key_default;
	_root_integer_4 _virtual_key_menu_0 = _virtual_key_default;
	_root_integer_4 _virtual_key_menu_1 = _virtual_key_default;
	_root_integer_4 _virtual_key_menu_2 = _virtual_key_default;
	_root_integer_4 _virtual_key_menu_3 = _virtual_key_default;
	_root_integer_4 _virtual_key_move_left = _virtual_key_default;
	_root_integer_4 _virtual_key_move_right = _virtual_key_default;
	_root_integer_4 _virtual_key_move_up = _virtual_key_default;
	_root_integer_4 _virtual_key_move_down = _virtual_key_default;
	_root_integer_4 _virtual_key_show_coordinates = _virtual_key_default;
	_root_integer_4_2 _position_default = (_root_integer_4_2){0, 0};
	_root_integer_4_2 _position_price = _position_default;
	_root_integer_4_2 _position_buy = _position_default;
	_root_integer_4_2 _position_buy_cancel = _position_default;
	_root_integer_4_2 _position_sell = _position_default;
	_root_integer_4_2 _position_sell_cancel = _position_default;
	_root_integer_4_2 _position_cancel = _position_default;
	_root_integer_4_2 _position_entry = _position_default;
	_root_integer_4_2 _position_exit = _position_default;
	_root_integer_4_2 _position_volume = _position_default;
	_root_integer_4_2 _position_volume_up = _position_default;
	_root_integer_4_2 _position_volume_down = _position_default;
	_root_integer_4_2 _position_scroll_up = _position_default;
	_root_integer_4_2 _position_scroll_down = _position_default;
	_root_integer_4 _menu_0_names_count = 0;
	_root_integer_1_unsigned *_menu_0_names = NULL;
	_root_integer_4 _menu_1_names_count = 0;
	_root_integer_1_unsigned *_menu_1_names = NULL;
	_root_integer_4 _menu_2_names_count = 0;
	_root_integer_1_unsigned *_menu_2_names = NULL;
	_root_integer_4 _menu_3_names_count = 0;
	_root_integer_1_unsigned *_menu_3_names = NULL;
	FILE *_file = fopen(u8"a.txt", u8"rb");
	if (!_file) {
		return 1;
	}
	_root_integer_4 _line_count = 1024;
	_root_character _line[_line_count];
	while (fgets(_line, _line_count, _file)) {
		_root_integer_4 _name_count = 1024;
		_root_character _name[_name_count];
		sscanf(_line, u8"%s", _name);
		if (!strcmp(_name, u8"_virtual_key_mode_0")) {
			sscanf(_line, u8"%*s%i", &_virtual_key_mode_0);
			continue;
		}
		if (!strcmp(_name, u8"_virtual_key_mode_1")) {
			sscanf(_line, u8"%*s%i", &_virtual_key_mode_1);
			continue;
		}
		if (!strcmp(_name, u8"_virtual_key_mode_2")) {
			sscanf(_line, u8"%*s%i", &_virtual_key_mode_2);
			continue;
		}
		if (!strcmp(_name, u8"_virtual_key_mode_3")) {
			sscanf(_line, u8"%*s%i", &_virtual_key_mode_3);
			continue;
		}
		if (!strcmp(_name, u8"_virtual_key_price_up")) {
			sscanf(_line, u8"%*s%i", &_virtual_key_price_up);
			continue;
		}
		if (!strcmp(_name, u8"_virtual_key_price_down")) {
			sscanf(_line, u8"%*s%i", &_virtual_key_price_down);
			continue;
		}
		if (!strcmp(_name, u8"_virtual_key_buy")) {
			sscanf(_line, u8"%*s%i", &_virtual_key_buy);
			continue;
		}
		if (!strcmp(_name, u8"_virtual_key_buy_cancel")) {
			sscanf(_line, u8"%*s%i", &_virtual_key_buy_cancel);
			continue;
		}
		if (!strcmp(_name, u8"_virtual_key_sell")) {
			sscanf(_line, u8"%*s%i", &_virtual_key_sell);
			continue;
		}
		if (!strcmp(_name, u8"_virtual_key_sell_cancel")) {
			sscanf(_line, u8"%*s%i", &_virtual_key_sell_cancel);
			continue;
		}
		if (!strcmp(_name, u8"_virtual_key_scroll_up")) {
			sscanf(_line, u8"%*s%i", &_virtual_key_scroll_up);
			continue;
		}
		if (!strcmp(_name, u8"_virtual_key_scroll_down")) {
			sscanf(_line, u8"%*s%i", &_virtual_key_scroll_down);
			continue;
		}
		if (!strcmp(_name, u8"_virtual_key_entry")) {
			sscanf(_line, u8"%*s%i", &_virtual_key_entry);
			continue;
		}
		if (!strcmp(_name, u8"_virtual_key_exit")) {
			sscanf(_line, u8"%*s%i", &_virtual_key_exit);
			continue;
		}
		if (!strcmp(_name, u8"_virtual_key_volume_up")) {
			sscanf(_line, u8"%*s%i", &_virtual_key_volume_up);
			continue;
		}
		if (!strcmp(_name, u8"_virtual_key_volume_down")) {
			sscanf(_line, u8"%*s%i", &_virtual_key_volume_down);
			continue;
		}
		if (!strcmp(_name, u8"_virtual_key_enter")) {
			sscanf(_line, u8"%*s%i", &_virtual_key_enter);
			continue;
		}
		if (!strcmp(_name, u8"_virtual_key_escape")) {
			sscanf(_line, u8"%*s%i", &_virtual_key_escape);
			continue;
		}
		if (!strcmp(_name, u8"_virtual_key_menu_0")) {
			sscanf(_line, u8"%*s%i", &_virtual_key_menu_0);
			continue;
		}
		if (!strcmp(_name, u8"_virtual_key_menu_1")) {
			sscanf(_line, u8"%*s%i", &_virtual_key_menu_1);
			continue;
		}
		if (!strcmp(_name, u8"_virtual_key_menu_2")) {
			sscanf(_line, u8"%*s%i", &_virtual_key_menu_2);
			continue;
		}
		if (!strcmp(_name, u8"_virtual_key_menu_3")) {
			sscanf(_line, u8"%*s%i", &_virtual_key_menu_3);
			continue;
		}
		if (!strcmp(_name, u8"_virtual_key_move_left")) {
			sscanf(_line, u8"%*s%i", &_virtual_key_move_left);
			continue;
		}
		if (!strcmp(_name, u8"_virtual_key_move_right")) {
			sscanf(_line, u8"%*s%i", &_virtual_key_move_right);
			continue;
		}
		if (!strcmp(_name, u8"_virtual_key_move_up")) {
			sscanf(_line, u8"%*s%i", &_virtual_key_move_up);
			continue;
		}
		if (!strcmp(_name, u8"_virtual_key_move_down")) {
			sscanf(_line, u8"%*s%i", &_virtual_key_move_down);
			continue;
		}
		if (!strcmp(_name, u8"_virtual_key_show_coordinates")) {
			sscanf(_line, u8"%*s%i", &_virtual_key_show_coordinates);
			continue;
		}
		if (!strcmp(_name, u8"_position_price")) {
			sscanf(_line, u8"%*s%i%i", &_position_price._0, &_position_price._1);
			continue;
		}
		if (!strcmp(_name, u8"_position_buy")) {
			sscanf(_line, u8"%*s%i%i", &_position_buy._0, &_position_buy._1);
			continue;
		}
		if (!strcmp(_name, u8"_position_buy_cancel")) {
			sscanf(_line, u8"%*s%i%i", &_position_buy_cancel._0, &_position_buy_cancel._1);
			continue;
		}
		if (!strcmp(_name, u8"_position_sell")) {
			sscanf(_line, u8"%*s%i%i", &_position_sell._0, &_position_sell._1);
			continue;
		}
		if (!strcmp(_name, u8"_position_sell_cancel")) {
			sscanf(_line, u8"%*s%i%i", &_position_sell_cancel._0, &_position_sell_cancel._1);
			continue;
		}
		if (!strcmp(_name, u8"_position_cancel")) {
			sscanf(_line, u8"%*s%i%i", &_position_cancel._0, &_position_cancel._1);
			continue;
		}
		if (!strcmp(_name, u8"_position_entry")) {
			sscanf(_line, u8"%*s%i%i", &_position_entry._0, &_position_entry._1);
			continue;
		}
		if (!strcmp(_name, u8"_position_exit")) {
			sscanf(_line, u8"%*s%i%i", &_position_exit._0, &_position_exit._1);
			continue;
		}
		if (!strcmp(_name, u8"_position_volume")) {
			sscanf(_line, u8"%*s%i%i", &_position_volume._0, &_position_volume._1);
			continue;
		}
		if (!strcmp(_name, u8"_position_volume_up")) {
			sscanf(_line, u8"%*s%i%i", &_position_volume_up._0, &_position_volume_up._1);
			continue;
		}
		if (!strcmp(_name, u8"_position_volume_down")) {
			sscanf(_line, u8"%*s%i%i", &_position_volume_down._0, &_position_volume_down._1);
			continue;
		}
		if (!strcmp(_name, u8"_position_scroll_up")) {
			sscanf(_line, u8"%*s%i%i", &_position_scroll_up._0, &_position_scroll_up._1);
			continue;
		}
		if (!strcmp(_name, u8"_position_scroll_down")) {
			sscanf(_line, u8"%*s%i%i", &_position_scroll_down._0, &_position_scroll_down._1);
			continue;
		}
		if (!strcmp(_name, u8"_menu_0_names")) {
			_root_integer_4 _index_0 = _root_string_count((_root_integer_1_unsigned *)u8"_menu_0_names") + 1;
			_root_integer_4 _index_1 = _root_string_find_last((_root_integer_1_unsigned *)_line, 0x09) + 1;
			_menu_0_names_count = _index_1 - _index_0;
			_menu_0_names = malloc(sizeof(_root_integer_1_unsigned) * _menu_0_names_count);
			memcpy(_menu_0_names, _line + _index_0, _menu_0_names_count);
			_root_string_replace(_menu_0_names_count, _menu_0_names, 0x09, 0x00);
			continue;
		}
		if (!strcmp(_name, u8"_menu_1_names")) {
			_root_integer_4 _index_0 = _root_string_count((_root_integer_1_unsigned *)u8"_menu_1_names") + 1;
			_root_integer_4 _index_1 = _root_string_find_last((_root_integer_1_unsigned *)_line, 0x09) + 1;
			_menu_1_names_count = _index_1 - _index_0;
			_menu_1_names = malloc(sizeof(_root_integer_1_unsigned) * _menu_1_names_count);
			memcpy(_menu_1_names, _line + _index_0, _menu_1_names_count);
			_root_string_replace(_menu_1_names_count, _menu_1_names, 0x09, 0x00);
			continue;
		}
		if (!strcmp(_name, u8"_menu_2_names")) {
			_root_integer_4 _index_0 = _root_string_count((_root_integer_1_unsigned *)u8"_menu_2_names") + 1;
			_root_integer_4 _index_1 = _root_string_find_last((_root_integer_1_unsigned *)_line, 0x09) + 1;
			_menu_2_names_count = _index_1 - _index_0;
			_menu_2_names = malloc(sizeof(_root_integer_1_unsigned) * _menu_2_names_count);
			memcpy(_menu_2_names, _line + _index_0, _menu_2_names_count);
			_root_string_replace(_menu_2_names_count, _menu_2_names, 0x09, 0x00);
			continue;
		}
		if (!strcmp(_name, u8"_menu_3_names")) {
			_root_integer_4 _index_0 = _root_string_count((_root_integer_1_unsigned *)u8"_menu_3_names") + 1;
			_root_integer_4 _index_1 = _root_string_find_last((_root_integer_1_unsigned *)_line, 0x09) + 1;
			_menu_3_names_count = _index_1 - _index_0;
			_menu_3_names = malloc(sizeof(_root_integer_1_unsigned) * _menu_3_names_count);
			memcpy(_menu_3_names, _line + _index_0, _menu_3_names_count);
			_root_string_replace(_menu_3_names_count, _menu_3_names, 0x09, 0x00);
			continue;
		}
	}
	fclose(_file);
	while (1) {
		if (_root_key_state(_virtual_key_mode_0)) {
			_mode = 0;
			Sleep(_delay);
			continue;
		}
		if (_root_key_state(_virtual_key_mode_1)) {
			_mode = 1;
			Sleep(_delay);
			continue;
		}
		switch (_mode) {
			case 0 : {
				Sleep(_delay);
				continue;
			} break;
			case 1 : {
				if (_root_key_state(_virtual_key_price_up)) {
					if (_row_index > 0) {
						_row_index += -1;
					}
					_root_cursor_position_set(_root_integer_4_2_add(_position_price, (_root_integer_4_2){0, _row_height * _row_index}));
					Sleep(_delay);
					continue;
				}
				if (_root_key_state(_virtual_key_price_down)) {
					if (_row_index < _row_count - 1) {
						_row_index += 1;
					}
					_root_cursor_position_set(_root_integer_4_2_add(_position_price, (_root_integer_4_2){0, _row_height * _row_index}));
					Sleep(_delay);
					continue;
				}
				if (_root_key_state(_virtual_key_buy)) {
					_root_cursor_position_set(_root_integer_4_2_add(_position_buy, (_root_integer_4_2){0, _row_height * _row_index}));
					_root_click(0);
					_root_cursor_position_set(_root_integer_4_2_add(_position_price, (_root_integer_4_2){0, _row_height * _row_index}));
					Sleep(_delay);
					continue;
				}
				if (_root_key_state(_virtual_key_buy_cancel)) {
					_root_cursor_position_set(_root_integer_4_2_add(_position_buy_cancel, (_root_integer_4_2){0, _row_height * _row_index}));
					_root_click(1);
					_root_cursor_position_set(_root_integer_4_2_add(_root_integer_4_2_add(_position_buy_cancel, (_root_integer_4_2){0, _row_height * _row_index}), _position_cancel));
					_root_click(0);
					_root_cursor_position_set(_root_integer_4_2_add(_position_price, (_root_integer_4_2){0, _row_height * _row_index}));
					Sleep(_delay);
					continue;
				}
				if (_root_key_state(_virtual_key_sell)) {
					_root_cursor_position_set(_root_integer_4_2_add(_position_sell, (_root_integer_4_2){0, _row_height * _row_index}));
					_root_click(0);
					_root_cursor_position_set(_root_integer_4_2_add(_position_price, (_root_integer_4_2){0, _row_height * _row_index}));
					Sleep(_delay);
					continue;
				}
				if (_root_key_state(_virtual_key_sell_cancel)) {
					_root_cursor_position_set(_root_integer_4_2_add(_position_sell_cancel, (_root_integer_4_2){0, _row_height * _row_index}));
					_root_click(1);
					_root_cursor_position_set(_root_integer_4_2_add(_root_integer_4_2_add(_position_sell_cancel, (_root_integer_4_2){0, _row_height * _row_index}), _position_cancel));
					_root_click(0);
					_root_cursor_position_set(_root_integer_4_2_add(_position_price, (_root_integer_4_2){0, _row_height * _row_index}));
					Sleep(_delay);
					continue;
				}
				if (_root_key_state(_virtual_key_scroll_up)) {
					_root_cursor_position_set(_position_scroll_up);
					_root_click(0);
					_root_cursor_position_set(_root_integer_4_2_add(_position_price, (_root_integer_4_2){0, _row_height * _row_index}));
					Sleep(_delay);
					continue;
				}
				if (_root_key_state(_virtual_key_scroll_down)) {
					_root_cursor_position_set(_position_scroll_down);
					_root_click(0);
					_root_cursor_position_set(_root_integer_4_2_add(_position_price, (_root_integer_4_2){0, _row_height * _row_index}));
					Sleep(_delay);
					continue;
				}
				if (_root_key_state(_virtual_key_entry)) {
					_root_cursor_position_set(_position_entry);
					_root_click(0);
					_root_cursor_position_set(_position_volume);
					_root_click(0);
					_root_click(0);
					_root_type(0x2e);
					_root_cursor_position_set(_root_integer_4_2_add(_position_price, (_root_integer_4_2){0, _row_height * _row_index}));
					Sleep(_delay);
					continue;
				}
				if (_root_key_state(_virtual_key_exit)) {
					_root_cursor_position_set(_position_exit);
					_root_click(0);
					_root_cursor_position_set(_position_volume);
					_root_click(0);
					_root_click(0);
					_root_type(0x2e);
					_root_cursor_position_set(_root_integer_4_2_add(_position_price, (_root_integer_4_2){0, _row_height * _row_index}));
					Sleep(_delay);
					continue;
				}
				if (_root_key_state(_virtual_key_volume_up)) {
					_root_cursor_position_set(_position_volume_up);
					_root_click(0);
					_root_cursor_position_set(_root_integer_4_2_add(_position_price, (_root_integer_4_2){0, _row_height * _row_index}));
					Sleep(_delay);
					continue;
				}
				if (_root_key_state(_virtual_key_volume_down)) {
					_root_cursor_position_set(_position_volume_down);
					_root_click(0);
					_root_cursor_position_set(_root_integer_4_2_add(_position_price, (_root_integer_4_2){0, _row_height * _row_index}));
					Sleep(_delay);
					continue;
				}
				if (_root_key_state(_virtual_key_enter)) {
					_root_type(0x0d);
					Sleep(_delay);
					continue;
				}
				if (_root_key_state(_virtual_key_escape)) {
					_root_type(0x1b);
					Sleep(_delay);
					continue;
				}
				if (_root_key_state(_virtual_key_menu_0)) {
					for (_root_integer_4 _index = 0; _index < _menu_0_names_count; _index += _root_string_count(_menu_0_names + _index) + 1) {
						_root_integer_4 _name_count = 256;
						_root_integer_2_unsigned _name[_name_count];
						_root_convert_utf16_utf8(_name, _name_count, _menu_0_names + _index);
						_root_window_activate(_name);
					}
					Sleep(_delay);
					continue;
				}
				if (_root_key_state(_virtual_key_menu_1)) {
					for (_root_integer_4 _index = 0; _index < _menu_1_names_count; _index += _root_string_count(_menu_1_names + _index) + 1) {
						_root_integer_4 _name_count = 256;
						_root_integer_2_unsigned _name[_name_count];
						_root_convert_utf16_utf8(_name, _name_count, _menu_1_names + _index);
						_root_window_activate(_name);
					}
					Sleep(_delay);
					continue;
				}
				if (_root_key_state(_virtual_key_menu_2)) {
					for (_root_integer_4 _index = 0; _index < _menu_2_names_count; _index += _root_string_count(_menu_2_names + _index) + 1) {
						_root_integer_4 _name_count = 256;
						_root_integer_2_unsigned _name[_name_count];
						_root_convert_utf16_utf8(_name, _name_count, _menu_2_names + _index);
						_root_window_activate(_name);
					}
					Sleep(_delay);
					continue;
				}
				if (_root_key_state(_virtual_key_menu_3)) {
					for (_root_integer_4 _index = 0; _index < _menu_3_names_count; _index += _root_string_count(_menu_3_names + _index) + 1) {
						_root_integer_4 _name_count = 256;
						_root_integer_2_unsigned _name[_name_count];
						_root_convert_utf16_utf8(_name, _name_count, _menu_3_names + _index);
						_root_window_activate(_name);
					}
					Sleep(_delay);
					continue;
				}
				if (_root_key_state(_virtual_key_move_left)) {
					_root_cursor_position_set(_root_integer_4_2_add(_root_cursor_position_get(), (_root_integer_4_2){-1, 0}));
					Sleep(_delay);
					continue;
				}
				if (_root_key_state(_virtual_key_move_right)) {
					_root_cursor_position_set(_root_integer_4_2_add(_root_cursor_position_get(), (_root_integer_4_2){1, 0}));
					Sleep(_delay);
					continue;
				}
				if (_root_key_state(_virtual_key_move_up)) {
					_root_cursor_position_set(_root_integer_4_2_add(_root_cursor_position_get(), (_root_integer_4_2){0, -1}));
					Sleep(_delay);
					continue;
				}
				if (_root_key_state(_virtual_key_move_down)) {
					_root_cursor_position_set(_root_integer_4_2_add(_root_cursor_position_get(), (_root_integer_4_2){0, 1}));
					Sleep(_delay);
					continue;
				}
				if (_root_key_state(_virtual_key_show_coordinates)) {
					_root_integer_4_2 _position = _root_cursor_position_get();
					printf(u8"%i %i \n", _position._0, _position._1);
					Sleep(_delay);
					continue;
				}
				Sleep(_delay);
				continue;
			}
			default : {
				Sleep(_delay);
				continue;
			}
		}
	}
	return 0;
}
_root_integer_4 main(_root_integer_4 _parameters_count, _root_character **_parameters) {
	return _root_main(_parameters_count, _parameters);
}
