#include <locale.h>
#include <iostream>
#include <cmath>
#include <ctime>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <complex>
#include <iomanip>
#include <algorithm>
#define PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062862

const int exponent_polynomial = 3; // старшая степень полиномов
const int GALOIS_FIELD = 3; // поле Галуа (простое число)
const int coefficient_for_generate = exponent_polynomial + 1; // количество коэффициентов полинома
int number_of_polynomials_0_and_1_exp = 0; // количество полиномов нулевой и первой степеней (по определению непреводимые)
int number_of_polynomials_0_exp = 0; // количество полиномов нулевой степени (чтобы не делить на них в дальнейшем, т.к. остаток от деления будет 0)
int counter_for_irreducible_polynomials = 0; // счетчик для строк массива, который хранит неприводимые полиномы

using namespace std;

int** array_allocation_int(int string, int columns) // ДИНАМИЧЕСКОЕ ВЫДЕЛЕНИЕ ПАМЯТИ ПОД ДВУМЕРНЫЙ МАССИВ int
{
	int** array = new int* [string];
	for (int i = 0; i < string; i++)
	{
		array[i] = new int[columns];
	}
	return array;
}

float** array_allocation_float(int string, int columns) // ДИНАМИЧЕСКОЕ ВЫДЕЛЕНИЕ ПАМЯТИ ПОД ДВУМЕРНЫЙ МАССИВ float
{
	float** array = new float* [string];
	for (int i = 0; i < string; i++)
	{
		array[i] = new float[columns];
	}
	return array;
}

double** array_allocation_double(int string, int columns) // ДИНАМИЧЕСКОЕ ВЫДЕЛЕНИЕ ПАМЯТИ ПОД ДВУМЕРНЫЙ МАССИВ double
{
	double** array = new double* [string];
	for (int i = 0; i < string; i++)
	{
		array[i] = new double[columns];
	}
	return array;
}

complex<double>** array_allocation_complex(int string, int columns) // ДИНАМИЧЕСКОЕ ВЫДЕЛЕНИЕ ПАМЯТИ ПОД ДВУМЕРНЫЙ МАССИВ complex
{
	complex <double>** array = new complex <double> *[string];
	for (int i = 0; i < string; i++)
	{
		array[i] = new complex <double>[columns];
	}
	return array;
}

void memory_cleaning(int** array, int string) // ОСВОБОЖДЕНИЕ ВЫДЕЛЕННОЙ ПАМЯТИ ПОД ДВУМЕРНЫЙ МАССИВ ТИПА int
{
	for (int i = 0; i < string; i++)
	{
		delete[] array[i];
	}
	delete[] array;
}

void memory_cleaning(float** array, int string) // ОСВОБОЖДЕНИЕ ВЫДЕЛЕННОЙ ПАМЯТИ ПОД ДВУМЕРНЫЙ МАССИВ ТИПА float
{
	for (int i = 0; i < string; i++)
	{
		delete[] array[i];
	}
	delete[] array;
}

void memory_cleaning(double** array, int string) // ОСВОБОЖДЕНИЕ ВЫДЕЛЕННОЙ ПАМЯТИ ПОД ДВУМЕРНЫЙ МАССИВ ТИПА double
{
	for (int i = 0; i < string; i++)
	{
		delete[] array[i];
	}
	delete[] array;
}

void memory_cleaning(complex<double>** array, int string) // ОСВОБОЖДЕНИЕ ВЫДЕЛЕННОЙ ПАМЯТИ ПОД ДВУМЕРНЫЙ МАССИВ ТИПА complex
{
	for (int i = 0; i < string; i++)
	{
		delete[] array[i];
	}
	delete[] array;
}

//int** array_of_all_polynomials = array_allocation_int(number_polynomials_not_zero, coefficient_for_generate); // массив для хранения всех возможных полиномов

int* left(int a[], int N) // СМЕЩЕНИЕ НА 1 ВЛЕВО (используется в функции деления полиномов)
{
	int temp = a[0];
	for (int i = 0; i < N - 1; i++)
	{
		a[i] = a[i + 1];
	}
	a[N - 1] = temp;
	return a;
}

void conclusion_polynom(vector<int>& polynom, int exp) // ВЫВОД ПОЛИНОМОВ В КОНСОЛЬ
{
	int counter = exp;
	for (int i = 0; i <= exp; i++)
	{
		if (polynom[i] != 0)
		{
			if (i == exp)
			{
				cout << polynom[i];
			}
			else
			{
				cout << polynom[i] << "x^" << counter << " + ";
			}
		}
		counter--;
	}
}

int exponent(int* polynom, int exp) // ОПРЕДЕЛЕНИЕ СТАРШИХ СТЕПЕНЕЙ ПОЛИНОМОВ
{
	int new_exp = exp;
	for (int i = 0; i <= exp; i++)
	{
		if (polynom[i] == 0)
			new_exp--;
		else
			break;
	}
	return new_exp;
}

int exponent(vector<int>& vec, int exp) // ОПРЕДЕЛЕНИЕ СТАРШИХ СТЕПЕНЕЙ ПОЛИНОМОВ ПРИ ПОМОЩИ ВЕКТОРА
{
	int new_exp = exp;
	for (int i = 0; i <= exp; i++)
	{
		if (vec[i] == 0)
			new_exp--;
		else
			break;
	}
	return new_exp;
}

int module_polynom(int polynom, int size, int mod) // КОЭФФИЦИЕНТЫ, ВЗЯТЫЕ ПО МОДУЛЮ (используется в функции деления полиномов)
{
	int result;
	if (polynom >= 0)
	{
		result = polynom % mod;
		return result;
	}
	else
	{
		while (polynom < 0)
		{
			polynom = polynom + mod;
		}
		result = polynom % mod;
		return result;
	}
}

void generate_coefficient_for_polynom(int coef, int mod, int number_polynomials, int** array_of_all_polynomials) // ГЕНЕРАЦИЯ ВСЕХ ВОЗМОЖНЫХ ПОЛИНОМОВ НАД ПОЛЕМ GF(3)
{
	int digits_left_to_generate = coef; // оставшаяся длина для генерации (сколько символов еще должно быть сгенерировано)
	static int coefficient_combination[coefficient_for_generate]; // для запоминания предыдущих состояний (префикса)
	static int top = 0; // уровень заполнения массива "coefficient_combination"
	static int step_for_generate = 0;
	if (digits_left_to_generate == 0) // сгенерированна одна из ветвей дерева вариантов (крайний случай)
	{
		int zero = 0;
		int positive = 0;
		for (int i = 0; i < top; i++)
		{
			if (coefficient_combination[i] == 0)
			{
				zero++;
			}
			else
			{
				positive++;
			}
		}
		if (zero != coefficient_for_generate) // проверка на 0, т.е. чтобы все коэффициенты не были равны 0 одновременно
		{
			for (int i = 0; i < top; i++)
			{
				array_of_all_polynomials[step_for_generate][i] = coefficient_combination[i];
			}
			step_for_generate++; // приращение счетчика для строк массива "array_of_all_polynomials"
		}
		if ((positive == 1) && (coefficient_combination[top - 1] != 0)) // полиномы нулевой степени
		{
			number_of_polynomials_0_exp++;
			number_of_polynomials_0_and_1_exp++;
		}
		if (((positive == 1) && (coefficient_combination[top - 2] != 0)) ||
			((positive == 2) && (coefficient_combination[top - 1] != 0) && (coefficient_combination[top - 2] != 0))) // полиномы первой степени
		{
			number_of_polynomials_0_and_1_exp++;
			counter_for_irreducible_polynomials++;
		}
	}
	else
	{
		for (int value = 0; value < GALOIS_FIELD; value++)
		{
			coefficient_combination[top++] = value; // дерево вариантов, начиная с value
			generate_coefficient_for_polynom(digits_left_to_generate - 1, mod, number_polynomials, array_of_all_polynomials);
			top--;
		}
	}
}

int division_polynom(int* first_polynom_from_main, vector<int>& second_polynom, int size, int exp, int step, int mod, int* irreducible_polynom, int step_devider, vector <vector<int> >& irreducible_polynom_vec) // ДЕЛЕНИЕ ПОЛИНОМОВ
{
	int break_division = 0; // индикатор приводимости полинома
	int* first_polynom = new int[size];
	for (int i = 0; i < size; i++)
	{
		first_polynom[i] = first_polynom_from_main[i];
	}
	int exp_first = exponent(first_polynom, exp);
	int exp_second = exponent(second_polynom, exp);
	int result = 0; // число, на которое умножается делитель перед вычитанием
	int exp_result = exp_first - exp_second;
	int divider;
	for (int i = 0; i < size; i++)
	{
		if (second_polynom[i] != 0) // нахождение старшего коэффициента делителя
		{
			divider = i;
			break;
		}
	}
	int new_exp_first = exp_first; // для приращения степени при выполнении каждой рекурсии
	int variable = 0; // значение, при умножении которого на старший коэффициент делителя, получится коэффициент старший коэффициент делимого (по модулю)
	if (first_polynom[step] != 0) // нахождение старшего коэффициента делимого
	{
		while (first_polynom[step] != variable) // нахождение необходимых коэффициентов для вычитания из делимого
		{
			result++;
			variable = (second_polynom[divider] * result);
			variable = variable % mod;
		}
	}
	else
	{
		int step_reserve = step;
		while (first_polynom[step_reserve] == 0 && new_exp_first > 0) // поиск первого (слева) не нулевого коэффициента делимого
		{
			step_reserve++;
		}
		if (first_polynom[step_reserve] != 0 && new_exp_first > 0)
		{
			while (first_polynom[step_reserve] != variable) // нахождение необходимых коэффициентов для вычитания из делимого
			{
				result++;
				variable = (second_polynom[divider] * result);
				variable = variable % mod;
			}
		}
		else
		{
			if (step_devider == 0)
			{
				counter_for_irreducible_polynomials++;
				for (int i = 0; i < 1; i++)
				{
					vector<int> row;
					for (int j = 0; j < size; j++)
					{
						row.push_back(irreducible_polynom[j]);
						//irreducible_polynomials[counter_for_irreducible_polynomials][j] = irreducible_polynom[j];
					}
					irreducible_polynom_vec.push_back(row);
				}
				return break_division;
			}
			return break_division;
		}
	}
	int value_shift = new_exp_first - exp_second; // на сколько нужно делать сдвиг влево
	const int SIZE_SUBTRAHEND = coefficient_for_generate;
	int subtrahend[SIZE_SUBTRAHEND];
	for (int i = 0; i < SIZE_SUBTRAHEND; i++)
	{
		subtrahend[i] = second_polynom[i] * result;
		subtrahend[i] = subtrahend[i] % mod;
	}
	int* shift = new int[SIZE_SUBTRAHEND]; // СДВИГ ВЛЕВО
	int sum = 0; // значение, которое определяет, не равен ли нулю остаток
	if (value_shift >= 0)
	{
		if (value_shift > 0)
		{
			for (int i = 0; i < value_shift; i++)
			{
				shift = left(subtrahend, SIZE_SUBTRAHEND);
			}
		}
		else
		{
			for (int i = 0; i < SIZE_SUBTRAHEND; i++)
			{
				shift[i] = subtrahend[i];
			}
		}
		for (int i = 0; i < SIZE_SUBTRAHEND; i++)
		{
			first_polynom[i] = (first_polynom[i] - shift[i]);
			first_polynom[i] = module_polynom(first_polynom[i], SIZE_SUBTRAHEND, mod); // нахождение коэффициентов по модулю
			if (first_polynom[i] == 0)
			{
				sum++;
			}
		}
	}
	if (sum == coefficient_for_generate) // проверка полинома на приводимость
	{
		break_division = 8; // индикатор приводимого полинома (если разделился на цело, то цикл for в main по j завершается)
		return break_division;
	}
	else if (sum != coefficient_for_generate && value_shift == 0) // проверка полинома на неприводимость 
	{
		if (step_devider == 0)
		{
			counter_for_irreducible_polynomials++;
			for (int i = 0; i < 1; i++)
			{
				vector<int> row;
				for (int j = 0; j < size; j++)
				{
					row.push_back(irreducible_polynom[j]);
				}
				irreducible_polynom_vec.push_back(row);
			}
			return break_division;
		}
	}
	exp_first--;
	step++;
	while (value_shift > 0)
	{
		if (exp_second > exp_first)
		{
			return break_division;
		}
		else
			return division_polynom(first_polynom, second_polynom, size, exp, step, mod, irreducible_polynom, step_devider, irreducible_polynom_vec);
	}
	delete[] first_polynom;
	delete[] shift;
}

void print_vector(int array[], int size) // ВЫВОД В КОНСОЛЬ ОДНОМЕРНЫХ МАССИВОВ
{
	for (int i = 0; i < size; i++)
	{
		cout << array[i];
	}
}

void print_matrix(int** array, int string, int column) // ВЫВОД В КОНСОЛЬ ДВУХМЕРНЫХ МАССИВОВ int
{
	for (int i = 0; i < string; i++)
	{
		for (int j = 0; j < column; j++)
		{
			cout << array[i][j];
		}
		cout << endl;
	}
}

void print_matrix(vector <vector<int>>& vec, int string, int column) // ВЫВОД В КОНСОЛЬ ДВУХМЕРНЫХ ВЕКТОРОВ
{
	for (int i = 0; i < string; i++)
	{
		for (int j = 0; j < column; j++)
		{
			cout << vec[i][j];
		}
		cout << endl;
	}
}

int* right(int array[], int size) // СМЕЩЕНИЕ НА 1 ЭЛЕМЕНТ ВПРАВО
{
	int temp = array[size - 1];
	for (int i = size - 1; i >= 0; i--)
	{
		array[i] = array[i - 1];
	}
	array[0] = temp;
	return array;
}

void right(vector <int>& array, int size) // СМЕЩЕНИЕ НА 1 ЭЛЕМЕНТ ВПРАВО
{
	int temp = array[size - 1];
	for (int i = size - 1; i >= 0; i--)
	{
		array[i] = array[i - 1];
	}
	array[0] = temp;
	//return array;
}

int** registers_addition(int array[], int size) // ОПРЕДЕЛЕНИЕ НЕНУЛЕВЫХ КОЭФФИЦИЕНТОВ МНОГОЧЛЕНА (КРОМЕ НУЛЕВОЙ СТРЕПЕНИ)
{
	int exp = size - 1; // чтобы не учитывать нулевую степень полинома
	int string = 2; // 2 строки: первая - номера ненулевых регистров, вторая - их коэффициенты
	int** draft_coefficients = array_allocation_int(string, exp); // т.к. количество ненулевых элементов неизвестно
	int size_no_zero = 0;
	for (int i = 0, j = exp; i < exp, j > 0; i++, j--) // два условия для корректного определения ненулевых коэффициентов (задом наперед)
	{
		if (array[i] != 0)
		{
			size_no_zero++;
			if (j != 0)
				draft_coefficients[0][i] = j - 1;
			else
				draft_coefficients[0][i] = j;
		}
		else
			draft_coefficients[0][i] = -1;
		draft_coefficients[1][i] = array[i];
	}
	int size_coefficients = size_no_zero + 1; // +1 чтобы вернуть из функции размерность массива, которая находится в первом элементе массива
	int** coefficients_and_registrers = array_allocation_int(string, size_coefficients);
	coefficients_and_registrers[0][0] = size_coefficients; // запись размерности массива, содержащего номера нужных регистров в первой строке
	coefficients_and_registrers[1][0] = string; // количество строк массива
	int step = 0;
	for (int i = 1; i < size_coefficients; i++) // создание массива, который хранит в своем первом элементе значение 
	{                                           // размерности массива, а в остальных - номера регистров, которые будут складываться
		for (int j = step; j < exp; j++)
		{
			if (draft_coefficients[0][j] >= 0)
			{
				coefficients_and_registrers[0][i] = draft_coefficients[0][j];
				coefficients_and_registrers[1][i] = draft_coefficients[1][j];
				step = j + 1;
				break;
			}
		}
	}
	return coefficients_and_registrers;
	memory_cleaning(coefficients_and_registrers, string);
	memory_cleaning(draft_coefficients, string);

}

int one_PSP_symbol(int shift_register[], int size_shift_register, int** register_numbers_AND_coefficients, int string, int column) // СИНТЕЗ ОДНОГО СИМВОЛА ПСП
{
	int sum = 0;
	for (int i = 0; i < size_shift_register; i++) // цикл проходит по каждому регистру
	{
		for (int j = 0; j < column; j++) // цикл проходит по каждому значению массива, в котором хранятся номера регистров для сложения
		{
			if (i == register_numbers_AND_coefficients[0][j]) // если текущий регистр должен участвовать в сложении
			{
				sum += (shift_register[i] * register_numbers_AND_coefficients[1][j]); // домножение на коэффициент ненулевого элемента полинома
				sum = sum % GALOIS_FIELD;
			}
		}
	}
	int one_symbol_PSP = shift_register[size_shift_register - 1];
	right(shift_register, size_shift_register);
	shift_register[0] = sum;
	return one_symbol_PSP;
}

int multiplicative_inverse(int value, int mod) // НАХОЖДЕНИЕ МУЛЬТИПЛИКАТИВНОЙ ОБРАТНОЙ ЧИСЛА
{
	int q, size = 0, i = 0;
	int z[1000];
	int modul = mod;
	while (value != 1 && mod != 1)
	{
		if (mod > value)
		{
			size++;
			q = mod / value;
			mod = mod - value * q;
			z[i] = q;
		}
		else
		{
			size++;
			q = value / mod;
			value = value - mod * q;
			z[i] = q;
		}
		i++;
	}
	int* y = new int[size];
	for (int i = 0; i < size; i++)
	{
		if (i == 0)
		{
			y[i] = 0 - 1 * z[i];
		}
		if (i == 1)
		{
			y[i] = 1 - (y[i - 1] * z[i]);
		}
		if (i != 0 && i != 1)
		{
			y[i] = y[i - 2] - (y[i - 1] * z[i]);
		}
	}
	while (y[size - 1] < 0)
	{
		y[size - 1] = y[size - 1] + modul;
	}
	y[size - 1] = y[size - 1] % modul;
	return y[size - 1];
	delete[] y;
}

int* generation_PSP(int polynom[], int size, int exp, int L) // СИНТЕЗ ОДНОЙ ПСП
{
	for (int i = 0; i < size; i++)
	{
		polynom[i] = polynom[i] * -1; // инверсия, чтобы избавиться от минуса в правой части уравнения
		while (polynom[i] < 0)
		{
			polynom[i] = polynom[i] + GALOIS_FIELD;
		}
		if (polynom[size - 1] != 1) // если коэффициент свободного члена отличен от единицы - умножить все
			// коэффициенты полинома на мультипликативное обратное коэффициента свободного члена
		{
			int mult_inverse = multiplicative_inverse(polynom[size - 1], GALOIS_FIELD); // определение мультипликативной
				// обратной коэффициента свободного члена
			polynom[i] = (polynom[i] * mult_inverse) % GALOIS_FIELD;
		}
	}
	int* shift_register = new int[exp]; // схема регистров
	shift_register[0] = 1;
	for (int i = 1; i < exp; i++)
	{
		shift_register[i] = 0;
	}
	int** register_numbers_begin = registers_addition(polynom, size);
	int string = register_numbers_begin[1][0];
	int size_nubers_register = register_numbers_begin[0][0]; // размерность массива, который хранит номера регистров для сложения (количество столбцов матрицы)
	int column = size_nubers_register - 1;  // новая размерность для массива, которы не будет хранить первый столбец
	int** register_numbers_AND_coefficients = array_allocation_int(string, column); // массив, хранящий номера регистров, которые будут складываться (не хранит размерность)
	for (int i = 0; i < string; i++)
	{
		for (int j = 0; j < column; j++)
		{
			register_numbers_AND_coefficients[i][j] = register_numbers_begin[i][j + 1];
		}
	}
	memory_cleaning(register_numbers_begin, string);
	int* PSP = new int[L]; // ПСП
	for (int i = 0; i < L; i++)
	{
		PSP[i] = one_PSP_symbol(shift_register, exp, register_numbers_AND_coefficients, string, column); // нужно передавать двухмерный массив в функцию
	}
	return PSP;
	memory_cleaning(register_numbers_AND_coefficients, string);
	delete[] shift_register;
	delete[] PSP;
}

void PSP_to_complex_vector(int PSP[], complex <double> complex_vector[], int duration_PSP) // ПРЕДСТАВЛЕНИЕ ПСП В КОМПЛЕКСНОМ ВИДЕ
{
	for (int i = 0; i < duration_PSP; i++)
	{
		complex <double> z(cos((2 * PI / GALOIS_FIELD) * PSP[i]), sin((2 * PI / GALOIS_FIELD) * PSP[i]));
		complex_vector[i] = z;
	}
}

double* crossСorrelation_function(int PSP_first[], int PSP_second[], int L) // ВЗАИМНОКОРРЕЛЯЦИОННАЯ ФУНКЦИЯ ДЛЯ ПСП
{
	complex <double>* complex_vector_first = new complex <double>[L];
	PSP_to_complex_vector(PSP_first, complex_vector_first, L);
	complex <double>* complex_vector_second = new complex <double>[L];
	PSP_to_complex_vector(PSP_second, complex_vector_second, L);
	complex <double>* crosscorrelation_function_complex = new complex <double>[L]; // значения ВКФ
	double* crosscorrelation_function = new double[L]; // вещественные части ВКФ
	for (int i = 0; i < L; i++)
	{
		crosscorrelation_function[i] = 0;
		crosscorrelation_function_complex[i] = (0, 0);
		for (int j = 0; j < L; j++)
		{
			crosscorrelation_function_complex[i] += real(complex_vector_first[j] * conj(complex_vector_second[(i + j) % L]));
		}
		crosscorrelation_function[i] = real(crosscorrelation_function_complex[i]);
	}
	return crosscorrelation_function;
	delete[] complex_vector_first;
	delete[] complex_vector_second;
	delete[] crosscorrelation_function_complex;
	delete[] crosscorrelation_function;
}

bool comparison(double first_value, double second_value, double accuracy) // функция сравнения двух чисел типа double
{
	bool check;
	double delta = abs(first_value - second_value);
	if (delta <= accuracy) // равны ли числа с требуемой точностью
		check = true;
	else
		check = false;
	return check;
}

bool determining_PSP_max_duration(double* ACF, int size_PSP) // ОПРЕДЕЛЕНИЕ ПСП МАКСИМАЛЬНОЙ ДЛИТЕЛЬНОСТИ НА ОСНОВЕ ЕЁ АКФ
{
	bool check = false; // индикатор ПСП с максимальным периодом
	int counter = 0; // счетчик "-1"
	double accuracy = 1e-9; // требуемая точность сравнения
	for (int i = 1; i < size_PSP; i++) // все элементы АКФ, начиная с 1 и до конца, должны быть равны -1 
	{
		bool equality = comparison(ACF[i], -1.0, accuracy); // для сравнения чисел типа double при помощи функции "comparison"
		if (equality)
		{
			counter++;
		}
	}
	if (counter == size_PSP - 1)
	{
		check = true;
	}
	return check;
}

double** cos_angle(int** good_PSP, int string, int duration_PSP)
{
	double** array_cos = array_allocation_double(string, duration_PSP); // массив, хранящий значиния 
	// косинусов углов между ПСП
	complex <double>** complex_matrix = array_allocation_complex(string, duration_PSP); // массив, 
	// хранящий ПСП в комплексной форме
	for (int i = 0; i < string; i++)
	{
		for (int j = 0; j < duration_PSP; j++)
		{
			complex <double> z(cos((2 * PI / GALOIS_FIELD) * good_PSP[i][j]), sin((2 * PI / GALOIS_FIELD) * good_PSP[i][j]));
			complex_matrix[i][j] = z;
		}
	}
	for (int i = 0; i < string; i++)
	{
		for (int j = 0; j < string; j++)
		{
			double cos_PSP = 0;
			complex <double> test = (0, 0);
			double numerator_sum = 0; // числитель формулы
			double mod_PSP_1 = 0;
			double mod_PSP_2 = 0;
			for (int k = 0; k < duration_PSP; k++)
			{
				test += complex_matrix[i][k] * conj(complex_matrix[j][k]);
				numerator_sum += real(complex_matrix[i][k] * conj(complex_matrix[j][k]));
				mod_PSP_1 += pow(abs(complex_matrix[i][k]), 2);
				mod_PSP_2 += pow(abs(complex_matrix[j][k]), 2);
			}
			mod_PSP_1 = sqrt(mod_PSP_1);
			mod_PSP_2 = sqrt(mod_PSP_2);
			array_cos[i][j] = numerator_sum / (mod_PSP_1 * mod_PSP_2);
		}
	}
	return array_cos;
	memory_cleaning(array_cos, string);
	memory_cleaning(complex_matrix, string);
}

void sorting(int** array, int string, int column) // СОРТИРОВКА СТОЛБЦОВ ПО УБЫВАНИЮ ЗНАЧЕНИЙ ПЕРВОЙ СТРОКИ
{
	for (int i = 0; i < column; i++)
	{
		for (int j = 0; j < column; j++)
		{
			if (array[0][i] > array[0][j])
			{
				swap(array[0][i], array[0][j]);
				swap(array[1][i], array[1][j]);
			}
		}
	}
}

int* selection(int** PSP, int quantity, int length, double** array_cos, double min_value_cos, double accuracy) // ВЫБОР КВАЗИОРТОГОНАЛЬНЫХ ПСП
{
	int string = 2; // количество строк массива содержащего кол-во минимумов и порядок ПСП
	int** quantity_min_value_cos_string = array_allocation_int(string, quantity); // количество минимумов каждой ПСП
	cout << endl << "Вектор кол-ва минимумов: ";
	for (int i = 0; i < quantity; i++) // вычисление кол-ва минимумов в каждой строке (для каждой ПСП)
	{
		quantity_min_value_cos_string[0][i] = 0;
		quantity_min_value_cos_string[1][i] = i; // запись номера ПСП
		for (int j = 0; j < quantity; j++)
		{
			bool equality = comparison(array_cos[i][j], min_value_cos, accuracy);
			if (equality)
			{
				quantity_min_value_cos_string[0][i]++;
			}
		}
		cout << quantity_min_value_cos_string[0][i];
	}
	double average_value_of_minimums = 0; // среднее значение кол-ва минимиумов
	int positions_max = 0; // ПСП, имеющая максимальное число минимумов
	int max_quantity_min = quantity_min_value_cos_string[0][positions_max]; // максимальное кол-во минимумов среди всех строк
	for (int i = 0; i < quantity; i++) // вычисление среднего значения кол-ва минимиумов
	{
		if (quantity_min_value_cos_string[0][i] > max_quantity_min)
		{
			max_quantity_min = quantity_min_value_cos_string[0][i];
			positions_max = i;
		}
		average_value_of_minimums += quantity_min_value_cos_string[0][i];
	}
	average_value_of_minimums = average_value_of_minimums / quantity;
	cout << endl << "Среднее значение количества минимумов: " << average_value_of_minimums;
	cout << endl << "Максимальное значение минимумов в одной строке: " << max_quantity_min;
	cout << endl << "ПСП, имеющая максимальное число минимумов: " << positions_max;
	sorting(quantity_min_value_cos_string, string, quantity);
	cout << endl << "Сортировка по убыванию (первая строка - кол-во минимумов; вторая - номер ПСП): " << endl;
	for (int i = 0; i < string; i++)
	{
		for (int j = 0; j < quantity; j++)
		{
			cout << quantity_min_value_cos_string[i][j];
		}
		cout << endl;
	}
	int number_potential_PSP = 0; // кол-во ПСП, кол-во минимумов которых не ниже среднего
	for (int i = 0; i < quantity; i++) // вычисление количества ПСП, кол-во минимумов которых не ниже среднего
	{
		if (quantity_min_value_cos_string[0][i] >= average_value_of_minimums)
		{
			number_potential_PSP++;
		}
	}
	cout << endl << "Кол-во ПСП, кол-во минимумов которых не ниже среднего: " << number_potential_PSP << " / " << quantity;
	int** positions_min_value_cos = array_allocation_int(quantity, quantity); // в пересечении с какой ПСП 
	// образуется минимум
	cout << endl << endl << "В пересечении с какой ПСП образуется минимум: " << endl;
	for (int i = 0; i < quantity; i++) // запись позиций минимумов в каждой строке (для каждой ПСП)
	{
		cout << "ПСП " << i << " образует минимумы в пересечении с ПСП: ";
		for (int j = 0; j < quantity; j++)
		{
			bool equality = comparison(array_cos[i][j], min_value_cos, accuracy);
			if (equality)
			{
				positions_min_value_cos[i][j] = j;
				cout << positions_min_value_cos[i][j] << " ";
			}
			else
			{
				positions_min_value_cos[i][j] = -1;
			}
		}
		cout << endl;
	}
	vector <int> list_quasi_orthogonal_PSP(1); // перечень квазиортогональных ПСП (с взаимными минимумами)
	list_quasi_orthogonal_PSP.at(0) = positions_max; // запись ПСП, которая обладает наибольшим числом минимумов
	cout << endl << "Из чего выбирать (сортировка по убыванию числа пересечений, образующих минимум): " << endl;
	for (int i = 0; i < quantity; i++)
	{
		int counter_for_list_quasi_orthogonal_PSP = 0; // счетчик номеров ПСП с взаимыными минимумами
		for (int j = 0; j < quantity; j++)
		{
			if (positions_min_value_cos[quantity_min_value_cos_string[1][i]][j] != -1) // проверка на наличие 
			// минимума в пересечении ПСП
			{
				cout << positions_min_value_cos[quantity_min_value_cos_string[1][i]][j] << " "; // упорядоченная 
				// запись позиций минимумов в каждой строке (для каждой ПСП по убыванию числа минимумов)
				for (int k = 0; k < list_quasi_orthogonal_PSP.size(); k++)
				{
					if (list_quasi_orthogonal_PSP[k] == positions_min_value_cos[quantity_min_value_cos_string[1][i]][j])
						// присутствуют ли все выписанные номера ПСП в вектор "list_quasi_orthogonal_PSP" в перечне 
						// минимумов номеров текущей ПСП (в строке массива "positions_min_value_cos")
					{
						counter_for_list_quasi_orthogonal_PSP++;
					}
					if (counter_for_list_quasi_orthogonal_PSP == list_quasi_orthogonal_PSP.size()) // если в
					// перечне минимумов текущей ПСП (в строке массива "positions_min_value_cos") имеются все, 
					// выписанные ПСП в вектор "list_quasi_orthogonal_PSP" номера ПСП (элементы данного вектора),
					// то номер строки "positions_min_value_cos" добавляется к "list_quasi_orthogonal_PSP"
					{
						list_quasi_orthogonal_PSP.push_back(quantity_min_value_cos_string[1][i]);
					}
				}
			}
		}
		cout << endl;
	}
	int* quasi_orthogonal_PSP = new int[list_quasi_orthogonal_PSP.size() + 1]; // массив, хранящий 
	// квазиортогональные ПСП
	quasi_orthogonal_PSP[0] = list_quasi_orthogonal_PSP.size(); // чтобы вернуть размерность массива в main
	cout << endl << "Номера ПСП с взаимными минимумами: ";
	for (int i = 0, j = 1; i < list_quasi_orthogonal_PSP.size(), j < list_quasi_orthogonal_PSP.size() + 1; i++, j++)
	{
		quasi_orthogonal_PSP[j] = list_quasi_orthogonal_PSP[i];
		cout << quasi_orthogonal_PSP[j];
	}
	return quasi_orthogonal_PSP;
	memory_cleaning(quantity_min_value_cos_string, string);
	memory_cleaning(positions_min_value_cos, quantity);
	delete[] quasi_orthogonal_PSP;
}

int factorial(int value)
{
	int fact = 1;
	for (int i = 0; i < value; i++)
	{
		fact += fact * i;
	}
	return fact;
}

int main()
{
	setlocale(LC_CTYPE, "rus");
	srand((unsigned)time(NULL));

	// ГЕНЕРИРОВАНИЕ НЕПРИВОДИМЫХ ПОЛИНОМОВ ЛЮБОЙ СТЕПЕНИ, НАД ЛЮБЫМ ПОЛЕМ
	vector <vector <int>> vector_irreducible(GALOIS_FIELD); // вектор, изначально содержащий все неприводимые полиномы 1-ой степени
	for (int i = 0; i < GALOIS_FIELD; i++)
	{
		vector_irreducible[i].resize(coefficient_for_generate);
	}
	int mod = GALOIS_FIELD;
	int exp = exponent_polynomial; // степень делимого полинома (и делителя, не достающие коэфициенты заполняются нулями слева)
	const int SIZE = coefficient_for_generate; // количество коэффициентов полинома
	int number_polynomials = pow(mod, coefficient_for_generate);  // количество вариантов полиномов (для функции генерации)
	int number_polynomials_not_zero = number_polynomials - 1; // количество вариантов полиномов без учета нулевого (для функции генерации) 
	int** array_of_all_polynomials = array_allocation_int(number_polynomials_not_zero, coefficient_for_generate); // массив для хранения всех возможных полиномов
	cout << "Количество возможных вариантов полиномов: " << number_polynomials_not_zero << endl;
	generate_coefficient_for_polynom(SIZE, mod, number_polynomials_not_zero, array_of_all_polynomials);
	for (int i = 0; i < number_polynomials_not_zero; i++)
	{
		for (int j = 0; j < SIZE; j++)
		{
			cout << array_of_all_polynomials[i][j] << " ";
		}
		cout << endl;
	}

	int** all_polynomials = array_allocation_int(number_polynomials_not_zero, coefficient_for_generate); // двумерный массив, содержащий все возможные полиномы
	// (для передачи в функцию деления)
	int** all_polynomials_for_check = array_allocation_int(number_polynomials_not_zero, coefficient_for_generate); // копия массива all_polynomials 
	//(для записи, внутри функции деления, в вектор неприводимых полиномов)	
	for (int i = 0; i < number_polynomials_not_zero; i++) // все возможные полиномы стпени exponent_polynomial
	{
		for (int j = 0; j < coefficient_for_generate; j++)
		{
			all_polynomials[i][j] = array_of_all_polynomials[i][j];
			all_polynomials_for_check[i][j] = all_polynomials[i][j];
		}
	}
	for (int i = number_of_polynomials_0_exp, k = 0; i < number_polynomials_not_zero, k < GALOIS_FIELD; i++, k++)
	{
		for (int j = 0; j < coefficient_for_generate; j++)
		{
			vector_irreducible[k][j] = all_polynomials[i][j];
		}
	}
	counter_for_irreducible_polynomials--;
	for (int i = number_of_polynomials_0_and_1_exp; i < number_polynomials_not_zero; i++) // передавать в функцию деления полиномы 2-ой и большей степени (что делить) 
	{
		int exp_first = exponent(all_polynomials[i], exp);
		int step_devider = vector_irreducible.size() - 1; // сколько раз полином не был разделен
		int break_division = 0;
		for (int j = 0; j < vector_irreducible.size(); j++)  // передавать в функцию деления полиномы 1-ой и большей степени (на что делить) 
		{
			if (break_division != 8) // если разделилось на цело хотя бы раз, то цикл прекращается
			{
				int exp_second = exponent(vector_irreducible[j], exp);
				int step = 0;
				break_division = division_polynom(all_polynomials[i], vector_irreducible[j], SIZE, exp, step, mod, all_polynomials_for_check[i], step_devider, vector_irreducible);
				step_devider--; // если цикл не обнулился
			}
			else
			{
				break;
			}
		}
	}
	memory_cleaning(all_polynomials, number_polynomials_not_zero);
	memory_cleaning(all_polynomials_for_check, number_polynomials_not_zero);
	cout << endl << "Список коэффициентов неприводимых полиномов:" << endl;
	print_matrix(vector_irreducible, vector_irreducible.size(), coefficient_for_generate);
	cout << endl << "Список неприводимых полиномов:" << endl;
	for (int i = 0; i < vector_irreducible.size(); i++)
	{
		conclusion_polynom(vector_irreducible[i], exp);
		cout << endl;
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// СИНТЕЗ ПСП
	int numbenumber_of_irreducible_polynomials = 0; // количество неприводимых полиномов степени exp
	for (int i = 0; i < vector_irreducible.size(); i++)
	{
		if (vector_irreducible[i][0] != 0)
		{
			numbenumber_of_irreducible_polynomials++;
		}
	}
	cout << endl << "Количество неприводимых полиномов " << exp << "-ой степени: " << numbenumber_of_irreducible_polynomials << endl;
	int begin_polynom_for_PSP = vector_irreducible.size() - numbenumber_of_irreducible_polynomials; // строка, откуда начинаются полиномы максимальной степени (exp)
	cout << endl << "Полиномы, на основе которых будут синтезироваться ПСП:" << endl;
	int** irreducible_polynomials_for_PSP = array_allocation_int(numbenumber_of_irreducible_polynomials, coefficient_for_generate); // массив, хранящий полиномы, на основе которых будут синтезироваться ПСП
	for (int i = 0, k = begin_polynom_for_PSP; i < vector_irreducible.size(), k < vector_irreducible.size(); i++, k++)
	{
		for (int j = 0; j < coefficient_for_generate; j++)
		{
			irreducible_polynomials_for_PSP[i][j] = vector_irreducible[k][j];
			cout << irreducible_polynomials_for_PSP[i][j] << " ";
		}
		cout << endl;
	}
	int duration_PSP = pow(GALOIS_FIELD, exp) - 1; // период ПСП на основе полиномов степени exp
	int** PSP = array_allocation_int(numbenumber_of_irreducible_polynomials, duration_PSP); // массив, хранящий ПСП
	for (int i = 0; i < numbenumber_of_irreducible_polynomials; i++)
	{
		PSP[i] = generation_PSP(irreducible_polynomials_for_PSP[i], coefficient_for_generate, exp, duration_PSP);
	}
	memory_cleaning(irreducible_polynomials_for_PSP, numbenumber_of_irreducible_polynomials);
	///////////////////////////////////////////////////////////////////////////////////////////////////////////

	// АКФ
	cout << endl;
	double** matrix_autocorrelation_function = array_allocation_double(numbenumber_of_irreducible_polynomials, duration_PSP); // массив, хранящий значения АКФ
	vector <vector <int>> good_PSP; // массив, хранящий ПСП с хорошими АКФ
	for (int i = 0; i < numbenumber_of_irreducible_polynomials; i++)
	{
		bool check_good_PSP = false;
		matrix_autocorrelation_function[i] = crossСorrelation_function(PSP[i], PSP[i], duration_PSP); // заполнение матрицы АКФ
		check_good_PSP = determining_PSP_max_duration(matrix_autocorrelation_function[i], duration_PSP);
		if (check_good_PSP)
		{
			vector <int> row;
			for (int j = 0; j < duration_PSP; j++)
			{
				row.push_back(PSP[i][j]);
			}
			good_PSP.push_back(row);
		}
	}
	int** good_PSP_array = array_allocation_int(good_PSP.size(), duration_PSP); // массив, содержащий ПСП макс. периода
	for (int i = 0; i < good_PSP.size(); i++)
	{
		for (int j = 0; j < duration_PSP; j++)
		{
			good_PSP_array[i][j] = good_PSP[i][j];
		}
	}
	cout << endl << endl << "Количество ПСП, имеющих лишь один пик АКФ: " << good_PSP.size() << endl;
	cout << endl << "ПСП, имеющие лишь один пик АКФ: " << endl;
	print_matrix(good_PSP_array, good_PSP.size(), duration_PSP);
	memory_cleaning(PSP, numbenumber_of_irreducible_polynomials);
	memory_cleaning(matrix_autocorrelation_function, numbenumber_of_irreducible_polynomials);
	///////////////////////////////////////////////////////////////////////////////////////////////////////

	// ЦИКЛИЧЕСКИЕ СДВИГИ
	double** array_cos = array_allocation_double(good_PSP.size(), good_PSP.size()); // массив, хранящий значиния косинусов углов между ПСП
	array_cos = cos_angle(good_PSP_array, good_PSP.size(), duration_PSP);
	int quantity_min_value_cos_begin = 0; // количество минимальных значений в исходной матрице
	double length_PSP_float = duration_PSP; // для получения требуемого результата при делении
	double min_value_cos = -1 / length_PSP_float;
	double accuracy = 1e-8; // требуемая точность сравнения значения матрицы косинусов с её минимальным значением
	cout << endl << "min_value_cos = " << min_value_cos << endl;
	cout << endl << "Матрица значений косиусов углов между ПСП:" << endl;
	for (int i = 0; i < good_PSP.size(); i++)
	{
		for (int j = 0; j < good_PSP.size(); j++)
		{
			cout << setw(10) << array_cos[i][j] << " ";
			bool equality = comparison(array_cos[i][j], min_value_cos, accuracy);
			if (equality)
			{
				quantity_min_value_cos_begin++;
			}
		}
		cout << endl;
	}
	cout << endl << "Количество минимальных элементов исходной матрицы: " << quantity_min_value_cos_begin << endl;
	int quantity_min_value_cos; // количество минимальных элементов после сдвига ПСП вправо на 1 элемент
	for (int i = 0; i < good_PSP.size(); i++)
	{
		for (int j = 0; j < good_PSP.size(); j++)
		{
			int position = 0; // сдвиг, при котором достигается максимальное количество элементов min_value_cos
			bool equality = comparison(array_cos[i][j], min_value_cos, accuracy);
			if (array_cos[i][j] != 1 && equality)
			{
				for (int k = 0; k < duration_PSP; k++)
				{
					quantity_min_value_cos = 0;
					good_PSP_array[j] = right(good_PSP_array[j], duration_PSP);
					array_cos = cos_angle(good_PSP_array, good_PSP.size(), duration_PSP);
					for (int n = 0; n < good_PSP.size(); n++)
					{
						for (int m = 0; m < good_PSP.size(); m++)
						{
							bool equality = comparison(array_cos[n][m], min_value_cos, accuracy);
							if (equality)
							{
								quantity_min_value_cos++;
							}
						}
					}
					if (quantity_min_value_cos > quantity_min_value_cos_begin)
					{
						quantity_min_value_cos_begin = quantity_min_value_cos;
						position = k;
					}
				}
				cout << endl << "Кол-во минимальных значений: " << quantity_min_value_cos_begin << endl;
				cout << "При позиции: " << position << endl;
				if (position != 0)
				{
					for (int k = 0; k < duration_PSP - 1 - position; k++)
					{
						good_PSP_array[j] = left(good_PSP_array[j], duration_PSP);
					}
				}
				for (int y = 0; y < duration_PSP; y++)
				{
					cout << good_PSP_array[j][y];
				}
				cout << endl;

				array_cos = cos_angle(good_PSP_array, good_PSP.size(), duration_PSP);
				cout << endl << "Матрица значений косиусов углов между ПСП:" << endl;
				for (int n = 0; n < good_PSP.size(); n++)
				{
					for (int m = 0; m < good_PSP.size(); m++)
					{
						cout << setw(10) << array_cos[n][m] << " ";
					}
					cout << endl;
				}
			}
		}
	}
	cout << endl << "Полученные ПСП: " << endl;
	print_matrix(good_PSP_array, good_PSP.size(), duration_PSP);
	//////////////////////////////////////////////////////////////////////////////////////////////////////

	// ВЫБОР КВАЗИОРТОГОНАЛЬНЫХ ПСП
	int* quasi_orthogonal_PSP = selection(good_PSP_array, good_PSP.size(), duration_PSP, array_cos, min_value_cos, accuracy);
	cout << endl << endl << "MAIN, quasi_orthogonal_PSP: ";
	int total_amount = quasi_orthogonal_PSP[0]; // сколько всего квазиортогональных ПСП
	for (int i = 1; i <= total_amount; i++)
	{
		cout << quasi_orthogonal_PSP[i];
	}
	int** quasi_orthogonal_PSP_full = array_allocation_int(total_amount, duration_PSP);
	cout << endl << endl << "Квазиортогональные ПСП:" << endl;
	for (int i = 1, k = 0; i <= total_amount, k < total_amount; i++, k++)
	{
		for (int j = 0; j < duration_PSP; j++)
		{
			quasi_orthogonal_PSP_full[k][j] = good_PSP_array[quasi_orthogonal_PSP[i]][j];
			cout << quasi_orthogonal_PSP_full[k][j];
		}
		cout << endl;
	}
	double** array_cos_quasi_orthogonal_PSP = array_allocation_double(total_amount, total_amount);
	array_cos_quasi_orthogonal_PSP = cos_angle(quasi_orthogonal_PSP_full, total_amount, duration_PSP);
	cout << endl << "Матрица значений косиусов углов между квазиортогональными ПСП:" << endl;
	for (int i = 0; i < total_amount; i++)
	{
		for (int j = 0; j < total_amount; j++)
		{
			cout << setw(10) << array_cos_quasi_orthogonal_PSP[i][j] << " ";
		}
		cout << endl;
	}
	//////////////////////////////////////////////////////////////////////////////////////////////////////

	// ВКФ КВАЗИОРТОГОНАЛЬНЫХ ПСП
	int quantity_PSP_for_CCF = 2; // сколько ПСП участвуют в ВКФ
	int quantity_combinations_CCF = factorial(total_amount) / (factorial(quantity_PSP_for_CCF) * factorial(total_amount - quantity_PSP_for_CCF)); // сколько всего неповторяющихся векторов ВКФ
	cout << endl << "Количество неповторяющихся векторов ВКФ: " << quantity_combinations_CCF << endl;
	double** matrix_cross_correlation_function = array_allocation_double(quantity_combinations_CCF, duration_PSP); // массив ВКФ
	cout << endl << "Элементы выше главной диагонали:" << endl;
	for (int i = 0; i < total_amount; i++)
	{
		for (int j = 0; j < total_amount; j++)
		{
			if (j > i)
				cout << setw(10) << array_cos_quasi_orthogonal_PSP[i][j] << " " << "(ПСП " << i << " с ПСП " << j << ")" << "\t";
		}
		cout << endl;
	}
	int counter_size_matrix_CCF = 0;
	for (int i = 0; i < total_amount; i++)
	{
		for (int j = 0; j < total_amount; j++)
		{
			if (j > i)
			{
				matrix_cross_correlation_function[counter_size_matrix_CCF] = crossСorrelation_function(quasi_orthogonal_PSP_full[i], quasi_orthogonal_PSP_full[j], duration_PSP); // заполнение матрицы ВКФ
				counter_size_matrix_CCF++;
			}
		}
	}
	cout << endl << "Результат вычисления ВКФ:" << endl;
	for (int i = 0; i < quantity_combinations_CCF; i++)
	{
		for (int j = 0; j < duration_PSP; j++)
		{
			cout << setw(4) << matrix_cross_correlation_function[i][j];
		}
		cout << endl;
	}
	double* max_CCF = new double[duration_PSP]; // вектор максимумов модулей ВКФ
	cout << endl << "Вектор максимумов модулей ВКФ:" << endl;
	for (int i = 0; i < duration_PSP; i++) // выбор максимумов модулей ВКФ
	{
		max_CCF[i] = abs(matrix_cross_correlation_function[0][i]);
		for (int j = 0; j < quantity_combinations_CCF; j++)
		{
			if (max_CCF[i] < abs(matrix_cross_correlation_function[j][i]))
			{
				max_CCF[i] = abs(matrix_cross_correlation_function[j][i]);
			}
		}
		cout << max_CCF[i] << " ";
	}
	double* min_CCF = new double[duration_PSP]; // вектор минимумов модулей ВКФ
	cout << endl << "Вектор минимумов модулей ВКФ:" << endl;
	for (int i = 0; i < duration_PSP; i++) // выбор максимумов ВКФ
	{
		min_CCF[i] = abs(matrix_cross_correlation_function[0][i]);
		for (int j = 0; j < quantity_combinations_CCF; j++)
		{
			if (min_CCF[i] > abs(matrix_cross_correlation_function[j][i]))
			{
				min_CCF[i] = abs(matrix_cross_correlation_function[j][i]);
			}
		}
		cout << min_CCF[i] << " ";
	}

	memory_cleaning(array_of_all_polynomials, number_polynomials_not_zero);
	memory_cleaning(good_PSP_array, good_PSP.size());
	memory_cleaning(array_cos, good_PSP.size());
	memory_cleaning(quasi_orthogonal_PSP_full, total_amount);
	memory_cleaning(array_cos_quasi_orthogonal_PSP, total_amount);
	memory_cleaning(matrix_cross_correlation_function, quantity_combinations_CCF);
	delete[] max_CCF;
	delete[] min_CCF;
	return 0;
}