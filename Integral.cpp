#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

const long double a = 1.1;
const long double b = 2.5;
const long double alfa = 0.4;
const long double betta = 0;
const long double pi = 3.14159265;
const long double eps = 0.000001;
const long double Integral1 = 14.27314090;
const long double Integral2 = 18.60294786;

long double f(long double x) {
	return 0.5 * cos(2 * x) * exp(0.4 * x) + 2.4 * sin(1.5 * x) * exp(-6 * x) + 6 * x;
}

long double p(long double x) {
	return pow((x - a), -alfa);
}

long double fp(long double x) {
	return f(x) * p(x);
}

long double middle(int n) {
	int i;
	long double c = (b - a) / n, x1, x2;
	long double F = 0;
	for (i = 0; i < n; i++) {
		x1 = a + i * c;
		x2 = x1 + c;
		F += (x2 - x1) * f((x1 + x2) / 2);
	}
	return F;
}

long double left(int n) {
	int i;
	long double c = (b - a) / n, x1, x2;
	long double F = 0;
	for (i = 0; i < n; i++) {
		x1 = a + i * c;
		x2 = x1 + c;
		F += (x2 - x1) * f(x1);
	}
	return F;
}

long double trapezoid(int n) {
	int i;
	long double c = (b - a) / n, x1, x2;
	long double F = 0;
	for (i = 0; i < n; i++) {
		x1 = a + i * c;
		x2 = x1 + c;
		F += (x2 - x1) * (f(x1) + f(x2)) / 2;
	}
	return F;
}

long double Simpson(int n) {
	int i;
	long double c = (b - a) / n, x1, x2;
	long double F = 0;
	for (i = 0; i < n; i++) {
		x1 = a + i * c;
		x2 = x1 + c;
		F += (x2 - x1) * (f(x1) + f(x2) + 4 * f((x1 + x2) / 2)) / 6;
	}
	return F;
}

long double Newton(int n) {
	int i;
	long double c = (b - a) / n, x1, x2, x12;
	long double nu0, nu1, nu2, A1, A2, A3;
	long double F = 0;
	for (i = 0; i < n; i++) {
		x1 = a + i * c;
		x2 = x1 + c;
		x12 = (x1 + x2) / 2;
		nu0 = (pow((x2 - a), 1 - alfa) - pow((x1 - a), 1 - alfa)) / (1 - alfa);
		nu1 = (pow((x2 - a), 2 - alfa) - pow((x1 - a), 2 - alfa)) / (2 - alfa) + a * nu0;
		nu2 = (pow((x2 - a), 3 - alfa) - pow((x1 - a), 3 - alfa)) / (3 - alfa) + 2 * a * nu1 - a * a * nu0;

		A1 = (nu2 - nu1 * (x12 + x2) + nu0 * x12 * x2) / ((x12 - x1) * (x2 - x1));
		A2 = -(nu2 - nu1 * (x1 + x2) + nu0 * x1 * x2) / ((x12 - x1) * (x2 - x12));
		A3 = (nu2 - nu1 * (x12 + x1) + nu0 * x12 * x1) / ((x2 - x12) * (x2 - x1));
		F += A1 * f(x1) + A2 * f(x12) + A3 * f(x2);
	}
	return F;
}

long double Gauss(int n) {
	int i, j, k, l;
	long double c = (b - a) / n, x1, x2, x12;
	long double nu0, nu1, nu2, nu3, nu4, nu5, p = 0, q = 0, r = 0, phi = 0, w;
	long double F = 0, d = 0;
	for (i = 0; i < n; i++) {
		x1 = a + i * c;
		x2 = x1 + c;
		nu0 = (pow((x2 - a), 1 - alfa) - pow((x1 - a), 1 - alfa)) / (1 - alfa);
		nu1 = (pow((x2 - a), 2 - alfa) - pow((x1 - a), 2 - alfa)) / (2 - alfa) + a * nu0;
		nu2 = (pow((x2 - a), 3 - alfa) - pow((x1 - a), 3 - alfa)) / (3 - alfa) + 2 * a * nu1 - a * a * nu0;
		nu3 = (pow((x2 - a), 4 - alfa) - pow((x1 - a), 4 - alfa)) / (4 - alfa) + 3 * a * nu2 - 3 * a * a * nu1 + pow(a, 3) * nu0;
		nu4 = (pow((x2 - a), 5 - alfa) - pow((x1 - a), 5 - alfa)) / (5 - alfa) + 4 * a * nu3 - 6 * a * a * nu2 + 4 * pow(a, 3) * nu1 - pow(a, 4) * nu0;
		nu5 = (pow((x2 - a), 6 - alfa) - pow((x1 - a), 6 - alfa)) / (6 - alfa) + 5 * a * nu4 - 10 * a * a * nu3 + 10 * pow(a, 3) * nu2 - 5 * pow(a, 4) * nu1 + pow(a, 5) * nu0;
		long double N[3][4] = { { nu0, nu1, nu2, -nu3 }, {nu1, nu2, nu3, -nu4 }, {nu2, nu3, nu4, -nu5} };
		long double aj[3] = { 0, 0, 0 };

		for (j = 0; j < 3; j++) {
			for (k = j + 1; k < 3; k++) {
				d = N[k][j] / N[j][j];
				for (l = j; l < 4; l++) {
					N[k][l] = N[k][l] - d * N[j][l];
				}
			}
		}

		aj[2] = N[2][3] / N[2][2];
		aj[1] = (N[1][3] - N[1][2] * aj[2]) / N[1][1];
		aj[0] = (N[0][3] - N[0][2] * aj[2] - N[0][1] * aj[1]) / N[0][0];

		q = (2 * pow(aj[2], 3)) / 54 - (aj[2] * aj[1]) / 6 + aj[0] / 2;
		p = (3 * aj[1] - aj[2] * aj[2]) / 9;

		if (q * q + p * p * p < 0) {
			r = sqrt(abs(p)) * abs(q) / q;
			phi = acos(q / pow(r, 3));
			w = aj[2];
			aj[0] = -2 * r * cos(phi / 3) - w / 3;
			aj[1] = 2 * r * cos(pi / 3 - phi / 3) - w / 3;
			aj[2] = 2 * r * cos(pi / 3 + phi / 3) - w / 3;
		}

		for (j = 0; j < 2; j++) {
			for (k = 0; k < 2; k++) {
				if (aj[k] > aj[k + 1]) {
					swap(aj[k], aj[k + 1]);
				}
			}
		}

		long double M[3][4] = { { 1, 1, 1, nu0 }, { aj[0], aj[1], aj[2], nu1 }, {aj[0] * aj[0], aj[1] * aj[1], aj[2] * aj[2], nu2} };
		long double Aj[3] = { 0, 0, 0 };

		for (j = 0; j < 3; j++) {
			for (k = j + 1; k < 3; k++) {
				d = M[k][j] / M[j][j];
				for (l = j; l < 4; l++) {
					M[k][l] = M[k][l] - d * M[j][l];
				}
			}
		}

		Aj[2] = M[2][3] / M[2][2];
		Aj[1] = (M[1][3] - M[1][2] * Aj[2]) / M[1][1];
		Aj[0] = (M[0][3] - M[0][2] * Aj[2] - M[0][1] * Aj[1]) / M[0][0];

		F += Aj[0] * f(aj[0]) + Aj[1] * f(aj[1]) + Aj[2] * f(aj[2]);
	}
	return F;
}

long double RichardsonN(int n) {
	int i, j, r = 0, k, l;
	long double Rh = 1, J, m, d, S;
	long double h = b - a;

	m = -log((Newton(4) - Newton(2)) / (Newton(2) - Newton(1))) / log(2);
	cout << m << endl;
	while (Rh > eps) {
		r++;
		long double* hi = new long double[r + 1];
		hi[0] = h;
		for (i = 1; i < r + 1; i++) {
			hi[i] = hi[i - 1] / 2;
		}
		long double** H = new long double* [r + 1];
		for (i = 0; i < r + 1; i++) {
			H[i] = new long double [r + 1];
		}
		for (i = 0; i < r + 1; i++) {
			H[i][0] = 1;
		}
		for (i = 0; i < r + 1; i++) {
			for (j = 1; j < r + 1; j++) {
				H[i][j] = -pow(hi[i], m + j - 1);
			}
		}
		long double* Sh = new long double[r + 1];
		for (i = 0; i < r + 1; i++) {
			Sh[i] = Newton(pow(2,i));
		}
		S = Sh[r];
		long double* C = new long double[r + 1];

		for (j = 0; j < r + 1; j++) {
			for (k = j + 1; k < r + 1; k++) {
				d = H[k][j] / H[j][j];
				for (l = j; l < r + 1; l++) {
					H[k][l] = H[k][l] - d * H[j][l];
				}
				Sh[k] = Sh[k] - d * Sh[j];
			}
		}

		for (j = r; j >= 0; j--) {
			for (k = j - 1; k >= 0; k--) {
				d = H[k][j] / H[j][j];
				for (l = k; l < r + 1; l++) {
					H[k][l] = H[k][l] - d * H[j][l];
				}
				Sh[k] = Sh[k] - d * Sh[j];
			}
		}

		for (i = 0; i < r + 1; i++) {
			C[i] = Sh[i] / H[i][i];
		}

		Rh = abs(C[0] - S);
		//cout << r << " " << Rh << endl;
	}
	cout << r << endl;
	h = h / pow(2, r);

	return h;
}

long double RichardsonG(int n) {
	int i, j, r = 0, k, l;
	long double Rh = 1, J, m, d, S;
	long double h = b - a;

	m = -log((Gauss(4) - Gauss(2)) / (Gauss(2) - Gauss(1))) / log(2);
	cout << m << endl;
	while (Rh > eps) {
		r++;
		long double* hi = new long double[r + 1];
		hi[0] = h;
		for (i = 1; i < r + 1; i++) {
			hi[i] = hi[i - 1] / 2;
		}
		long double** H = new long double* [r + 1];
		for (i = 0; i < r + 1; i++) {
			H[i] = new long double[r + 1];
		}
		for (i = 0; i < r + 1; i++) {
			H[i][0] = 1;
		}
		for (i = 0; i < r + 1; i++) {
			for (j = 1; j < r + 1; j++) {
				H[i][j] = -pow(hi[i], m + j - 1);
			}
		}
		long double* Sh = new long double[r + 1];
		for (i = 0; i < r + 1; i++) {
			Sh[i] = Gauss(pow(2, i));
		}
		S = Sh[r];
		long double* C = new long double[r + 1];

		for (j = 0; j < r + 1; j++) {
			for (k = j + 1; k < r + 1; k++) {
				d = H[k][j] / H[j][j];
				for (l = j; l < r + 1; l++) {
					H[k][l] = H[k][l] - d * H[j][l];
				}
				Sh[k] = Sh[k] - d * Sh[j];
			}
		}

		for (j = r; j >= 0; j--) {
			for (k = j - 1; k >= 0; k--) {
				d = H[k][j] / H[j][j];
				for (l = k; l < r + 1; l++) {
					H[k][l] = H[k][l] - d * H[j][l];
				}
				Sh[k] = Sh[k] - d * Sh[j];
			}
		}

		for (i = 0; i < r + 1; i++) {
			C[i] = Sh[i] / H[i][i];
		}

		Rh = abs(C[0] - S);
		//cout << r << " " << Rh << endl;
	}
	cout << r << endl;

	h = h / pow(2, r);

	return h;
}

long double RichardsonOpt(int n) {
	int i, j, r = 0, k, l;
	long double Rh = 1, J, m, d, S;
	long double h = b - a;
	long double h1 = h, h2 = h / 2, Rh1, hopt;
	m = -log((Gauss(4) - Gauss(2)) / (Gauss(2) - Gauss(1))) / log(2);
	cout << m << endl;
	Rh1 = (Gauss(2) - Gauss(1)) / (1 - pow(2, m));
	hopt = h1 * pow((eps / abs(Rh1)), 1 / m);
	h = hopt;
	cout << "hopt: " << hopt << endl;
	while (Rh > eps) {
		r++;
		long double* hi = new long double[r + 1];
		hi[0] = h;
		for (i = 1; i < r + 1; i++) {
			hi[i] = hi[i - 1] / 2;
		}
		long double** H = new long double* [r + 1];
		for (i = 0; i < r + 1; i++) {
			H[i] = new long double[r + 1];
		}
		for (i = 0; i < r + 1; i++) {
			H[i][0] = 1;
		}
		for (i = 0; i < r + 1; i++) {
			for (j = 1; j < r + 1; j++) {
				H[i][j] = -pow(hi[i], m + j - 1);
			}
		}
		long double* Sh = new long double[r + 1];
		for (i = 0; i < r + 1; i++) {
			Sh[i] = Gauss(pow(2, i));
		}
		S = Sh[r];
		long double* C = new long double[r + 1];

		for (j = 0; j < r + 1; j++) {
			for (k = j + 1; k < r + 1; k++) {
				d = H[k][j] / H[j][j];
				for (l = j; l < r + 1; l++) {
					H[k][l] = H[k][l] - d * H[j][l];
				}
				Sh[k] = Sh[k] - d * Sh[j];
			}
		}

		for (j = r; j >= 0; j--) {
			for (k = j - 1; k >= 0; k--) {
				d = H[k][j] / H[j][j];
				for (l = k; l < r + 1; l++) {
					H[k][l] = H[k][l] - d * H[j][l];
				}
				Sh[k] = Sh[k] - d * Sh[j];
			}
		}

		for (i = 0; i < r + 1; i++) {
			C[i] = Sh[i] / H[i][i];
		}

		Rh = abs(C[0] - S);
	}

	h = h / pow(2, r);

	return h;
}

int main() {
	setlocale(LC_ALL, "Rus");
	int i, j, k, q, l;

	ofstream fout1, fout2, fout3, fout4, fout5, fout6;
	fout1.open("text1.txt");
	fout2.open("text2.txt");
	fout3.open("text3.txt");
	fout4.open("text4.txt");
	fout5.open("text5.txt");
	fout6.open("text6.txt");

	cout << "Выберите часть задания" << endl;
	cout << "1. Квадратурные формулы Ньютона-Кот(е)са и Гаусса" << endl;
	cout << "2. Методы оценки погрешности составных квадратурных формул" << endl;
	cin >> q;
	if (q == 1) {
		cout << "Какие использовать составные квадратурные формулы?" << endl;
		cout << "1. Квадратурная формула среднего прямоугольника" << endl;
		cout << "2. Квадратурная формула левого прямоугольника" << endl;
		cout << "3. Квадратурная формула трапеции" << endl;
		cout << "4. Квадратурная формула Симпсона" << endl;
		cout << "5. Квадратурная формула Ньютона-Котса" << endl;
		cout << "6. Квадратурная формула Гаусса" << endl;
 
		cin >> l;
		
		switch (l) {
		case 1:
			for (i = 1; i < 100; i++) {
				cout << i << ": ";
				printf("%.8lf\n", middle(i));
				fout1 << i << " " << 10 * abs(middle(i) - Integral1) << endl;
			}
			break;
		case 2:
			for (i = 1; i < 100; i++) {
				cout << i << ": ";
				printf("%.8lf\n", left(i));
				fout2 << i << " " << 10 * abs(left(i) - Integral1) << endl;
			}
			break;
		case 3:
			for (i = 1; i < 100; i++) {
				cout << i << ": ";
				printf("%.8lf\n", trapezoid(i));
				fout3 << i << " " << 10 * abs(trapezoid(i) - Integral1) << endl;
			}
			break;
		case 4:
			for (i = 1; i < 100; i++) {
				cout << i << ": ";
				printf("%.8lf\n", Simpson(i));
				fout4 << i << " " << 1000000 * abs(Simpson(i) - Integral1) << endl;
			}
			break;
		case 5:
			for (i = 1; i < 200; i++) {
				cout << i << ": ";
				printf("%.10lf\n", Newton(i));
				fout5 << i << " " << 100000 * abs(Newton(i) - Integral2) << endl;
			}
			break;

		case 6:
			for (i = 1; i < 100; i++) {
				cout << i << ": ";
				printf("%.8lf\n", Gauss(i));
				fout6 << i << " " << 100000 * abs(Gauss(i) - Integral2) << endl;
			}
			break;
		}
	}
	if (q == 2) {
		cout << "Оценка погрешности методом Ричардсона" << endl;
		cout << "1. Для квадратурной формулы Ньютона-Котса" << endl;
		cout << "2. Для квадратурной формулы Гаусса" << endl;
		cout << "3. Для квадратурной формулы Гаусса и оптимального шага" << endl;
		cin >> l;
		switch (l) {
		case 1:
			cout << RichardsonN(10);
			break;
		case 2:
			cout << RichardsonG(10);
			break;
		case 3:
			cout << RichardsonOpt(10);
			break;
		}
	}

	return 0;
}