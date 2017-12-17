#include "stdafx.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace std;

class gFEM3D8Node {

private:
	// Он держит меня в подвале!!1 ПОМОГИТЕ!
	
	int n; //Êîëè÷åñòâî óçëîâ
	int p; //Êîëè÷åñòâî ñòåïåíåé ñâîáîäû
	int q, qg;

	int i, j, k, s;

	//<Ìàòåðèàëû
	double mu, lam;
	//Ìàòåðèàëû>

	//<Ìàòðèöû æåñòêîñòè
	Eigen::MatrixXd KL;							//ËÌÆ						//Êàê óäàëÿòü?
	Eigen::SparseMatrix<double> SparseKG;		//ÃÌÆ
												//Ìàòðèöû æåñòêîñòè>

												//<Ñåòêà
	double** Nodes;				//Êîîðäèíàòû óçëîâ (ãëîáàëüíûå)
	int** Elems;				//Ãëîáàëüíàÿ íóìåðàöèÿ óçëîâ ÊÝ 
	int NNodes;					//Êîë-âî óçëîâ
	int NElems;					//Êîë-âî ýëåìåíòîâ
	Eigen::VectorXd BCForce;	//Íàãðóçêè
	std::vector<int> BCKinem;	//Çàêðåïëåíèÿ (íîìåðà çàêðåïëåííûõ óçëîâ)
								//Ñåòêà>

								//<Ïàðàìåòðû äëÿ CT
	int m1;						//Êîë-âî âîêñåëåé ïî îñè x
	int m2;						//Êîë-âî âîêñåëåé ïî îñè y
	int m3;						//Êîë-âî âîêñåëåé ïî îñè z
	int h;						//Ñìåùåíèå â áèíàðíîì ôàéëå "ïî ÊÝ" "âíèç"
	double voxX, voxY, voxZ;	//Ðàçìåðû âîêñåëÿ
	double a;					//Ôèç. ðàçìåð ÊÝ ïî îñè x
	double b;					//Ôèç. ðàçìåð ÊÝ ïî îñè y
	double c;					//Ôèç. ðàçìåð ÊÝ ïî îñè z
	int** CTData;				//Áèíàðíûé ìàññèâ
								//Ïàðàìåòðû äëÿ CT>

								//<Ïåðåìåùåíèÿ
	Eigen::VectorXd Delta;
	//Ïåðåìåùåíèÿ>

public:

	gFEM3D8Node(void) {

		n = 8;
		p = 3;
		q = n*p;

		this->ReadMain();

		KL.resize(q, q);

	}

	~gFEM3D8Node(void)
	{

		for (i = 0; i < h; i++)
			delete CTData[i];

		delete CTData;

		//KL.operator delete;
		//KL.operator delete[];

	}

private:
	//Ñ÷èòûâàíèå îñíîâíîãî ôàéëà:
	void ReadMain() {

		//<Ïàðàìåòðû ìàòåðèàëà
		double EE = 200000; //[E] = ÌÏà, ò.ê.[vox] = ìì.;
		double nu = 0.28;
		mu = EE / (1 + nu) / 2;
		lam = 2 * mu * nu / (1 - 2 * nu);
		//Ïàðàìåòðû ìàòåðèàëà>

		NElems = 2; //Êîëè÷åñòâî ýëåìåíòîâ
		NNodes = 12; //Êîëè÷åñòâî óçëîâ

		qg = NNodes*p;

		BCForce.resize(qg);
		Delta.resize(qg);

		//<Ãðàíè÷íûå óñëîâèÿ:
		double F = 1200, sign = -1;

		//<Íàãðóçêè
		/*Ïðèêëàäûâàåì -Py ê 3,4,5,9,10 è 11 óçëàì.*/
		for (i = 0; i < qg; i++)
			BCForce(i) = 0;

		for (i = 3; i <= 5; i++)
			for (j = 0; j <= 6; j += 6)
				BCForce(i + j) = sign*F / 6;
		//Íàãðóçêè>

		std::cout << BCForce << std::endl;
		cout << 2017 << endl;

		//<Çàêðåïëåíèÿ
		/*Çàêðïåëÿåì 0,1,2,6,7 è 8 óçëû ALLDOF.*/
		for (i = 0; i <= 6; i += 6)
			for (j = 0; j <= 2; j++)
				for (k = 0; k <= 2 * NNodes; k += NNodes) {

					BCForce(i + j + k) = 0;
					BCKinem.push_back(i + j + k);

				}
		//Çàêðåïëåíèÿ>

		Nodes = new double*[NNodes];
		for (i = 0; i < NNodes; i++)
			Nodes[i] = new double[p];

		Elems = new int*[NElems];
		for (i = 0; i < NElems; i++)
			Elems[i] = new int[n];

		this->ReadBin();

		//<Ôèç. ðàçìåðíîñòü îáúåìà
		a = m1*voxX,
			b = m2*voxY,
			c = m3*voxZ;
		//Ôèç. ðàçìåðíîñòü îáúåìà>

		//<Êîîðäèíàòû óçëîâ
		/*Nodes[a1][b1]:
		a1 - íîìåð óçëà (ãëîáàëüíûé),
		a2 - íîìåð êîîðäèíàòû.*/

		double r = 0.;

		for (i = 0; i < NNodes / 4; i++) {

			Nodes[i][0] = r;
			Nodes[i][1] = 0;
			Nodes[i][2] = 0;

			Nodes[i + NNodes / 4][0] = r;
			Nodes[i + NNodes / 4][1] = b;
			Nodes[i + NNodes / 4][2] = 0;

			Nodes[i + NNodes / 2][0] = r;
			Nodes[i + NNodes / 2][1] = 0;
			Nodes[i + NNodes / 2][2] = c;

			Nodes[i + 3 * NNodes / 4][0] = r;
			Nodes[i + 3 * NNodes / 4][1] = b;
			Nodes[i + 3 * NNodes / 4][2] = c;

			r += a;

		}
		//Êîîðäèíàòû óçëîâ>

		//<Íîìåðà ýëåìåíòîâ â ÊÝ
		for (i = 0; i < 2; i++) {

			Elems[i][0] = i + 0;
			Elems[i][1] = i + 1;
			Elems[i][2] = i + 3;
			Elems[i][3] = i + 4;
			Elems[i][4] = i + 6;
			Elems[i][5] = i + 7;
			Elems[i][6] = i + 9;
			Elems[i][7] = i + 10;

		}
		//Íîìåðà ýëåìåíòîâ â ÊÝ>

	}

	//Ñ÷èòûâàíèå áèíàðíîãî ìàññèâà:
	void ReadBin() {

		//Ðàçìåðíîñòü âîêñåëÿ:
		voxX = 0.2;
		voxY = 0.2;
		voxZ = 0.2;

		//Ðàçìåðíîñòü CT:
		m1 = 5;
		m2 = 5;
		m3 = 5;
		h = m2*m3;

		//ÒÐÅÁÓÅÒÑß ÈÑÊËÞ×ÈÒÜ ÎÏÐÅÄÅËÅÍÈÅ ÐÀÇÌÅÐÍÎÑÒÈ ÌÀÑÑÈÂÀ! Â ÝÒÎÌ ÍÅÎÁÕÎÄÈÌÎÑÒÈ ÁÎËÜØÅ ÍÅÒ,
		//Ò.Ê. ÌÛ ÈÕ Ñ×ÈÒÛÂÀÅÌ ÈÇ ÔÀÉËÀ!

		//<Ñ÷èòûâàíèå â ìàññèâ áèíàðèçîâàííîé ìàòðèöû
		ifstream in("C:/Users/imm/Documents/Visual Studio 2017/Projects/ConsoleApplication1/BINC.txt");

		if (in.is_open()) {

			int count = 0;
			int temp;

			while (!in.eof()) {
				in >> temp;
				count++;
			}

			in.seekg(0, ios::beg);
			in.clear();

			int count_space = 0;
			char symbol;
			while (!in.eof()) {
				in.get(symbol);
				if (symbol == ' ') count_space++;
				if (symbol == '\n') break;
			}

			in.seekg(0, ios::beg);
			in.clear();

			//int s1 = count / (count_space + 1);
			//int s2 = count_space + 1;

			int s1 = h*NElems; //Ñòðîêè
			int s2 = m1; //Ñòîëáöû

			CTData = new int*[s1];
			for (i = 0; i < s1; i++)
				CTData[i] = new int[s2];

			for (i = 0; i < s1; i++)
				for (j = 0; j < s2; j++)
					in >> CTData[i][j];

			in.close();
		}

		else
			cout << "Ôàéë íå íàéäåí.";
		//Ñ÷èòûâàíèå â ìàññèâ áèíàðèçîâàííîé ìàòðèöû>

	}

	//Ïîñòðîåíèå ËÌÆ:
	void LocalMatrix() {

		Eigen::MatrixXd SYS, K, A, AT;

		SYS.resize(q, q);
		K.resize(q, q);
		A.resize(q, q);
		AT.resize(q, q);

		//<Èíòåãðèðîâàíèå ïî áèíàðèçîâàííîìó ìàññèâó
		int u1, u2;
		double x, y, z;

		for (u1 = 0; u1 < q; u1++)
			for (u2 = 0; u2 < q; u2++)
				K(u1, u2) = 0;

		for (k = 0; k < m3; k++)
			for (j = 0; j < m2; j++)
				for (i = 0; i < m1; i++) {

					x = (i + 0.5)*voxX;
					y = (j + 0.5)*voxY;
					z = (k + 0.5)*voxZ;

					//<Ïåðåîïðåäåëåíèå SYS
					SYS(0, 0) = 0;
					SYS(0, 1) = 0;
					SYS(0, 2) = 0;
					SYS(0, 3) = 0;
					SYS(0, 4) = 0;
					SYS(0, 5) = 0;
					SYS(0, 6) = 0;
					SYS(0, 7) = 0;
					SYS(0, 8) = 0;
					SYS(0, 9) = 0;
					SYS(0, 10) = 0;
					SYS(0, 11) = 0;
					SYS(0, 12) = 0;
					SYS(0, 13) = 0;
					SYS(0, 14) = 0;
					SYS(0, 15) = 0;
					SYS(0, 16) = 0;
					SYS(0, 17) = 0;
					SYS(0, 18) = 0;
					SYS(0, 19) = 0;
					SYS(0, 20) = 0;
					SYS(0, 21) = 0;
					SYS(0, 22) = 0;
					SYS(0, 23) = 0;
					SYS(1, 0) = 0;
					SYS(1, 1) = lam + 2 * mu;
					SYS(1, 2) = 0;
					SYS(1, 3) = 0;
					SYS(1, 4) = y*(lam + 2 * mu);
					SYS(1, 5) = 0;
					SYS(1, 6) = z*(lam + 2 * mu);
					SYS(1, 7) = y*z*(lam + 2 * mu);
					SYS(1, 8) = 0;
					SYS(1, 9) = 0;
					SYS(1, 10) = lam;
					SYS(1, 11) = 0;
					SYS(1, 12) = lam*x;
					SYS(1, 13) = lam*z;
					SYS(1, 14) = 0;
					SYS(1, 15) = lam*x*z;
					SYS(1, 16) = 0;
					SYS(1, 17) = 0;
					SYS(1, 18) = 0;
					SYS(1, 19) = lam;
					SYS(1, 20) = 0;
					SYS(1, 21) = lam*y;
					SYS(1, 22) = lam*x;
					SYS(1, 23) = lam*x*y;
					SYS(2, 0) = 0;
					SYS(2, 1) = 0;
					SYS(2, 2) = mu;
					SYS(2, 3) = 0;
					SYS(2, 4) = mu*x;
					SYS(2, 5) = mu*z;
					SYS(2, 6) = 0;
					SYS(2, 7) = mu*x*z;
					SYS(2, 8) = 0;
					SYS(2, 9) = mu;
					SYS(2, 10) = 0;
					SYS(2, 11) = 0;
					SYS(2, 12) = mu*y;
					SYS(2, 13) = 0;
					SYS(2, 14) = mu*z;
					SYS(2, 15) = mu*y*z;
					SYS(2, 16) = 0;
					SYS(2, 17) = 0;
					SYS(2, 18) = 0;
					SYS(2, 19) = 0;
					SYS(2, 20) = 0;
					SYS(2, 21) = 0;
					SYS(2, 22) = 0;
					SYS(2, 23) = 0;
					SYS(3, 0) = 0;
					SYS(3, 1) = 0;
					SYS(3, 2) = 0;
					SYS(3, 3) = mu;
					SYS(3, 4) = 0;
					SYS(3, 5) = mu*y;
					SYS(3, 6) = mu*x;
					SYS(3, 7) = mu*x*y;
					SYS(3, 8) = 0;
					SYS(3, 9) = 0;
					SYS(3, 10) = 0;
					SYS(3, 11) = 0;
					SYS(3, 12) = 0;
					SYS(3, 13) = 0;
					SYS(3, 14) = 0;
					SYS(3, 15) = 0;
					SYS(3, 16) = 0;
					SYS(3, 17) = mu;
					SYS(3, 18) = 0;
					SYS(3, 19) = 0;
					SYS(3, 20) = mu*y;
					SYS(3, 21) = 0;
					SYS(3, 22) = mu*z;
					SYS(3, 23) = mu*y*z;
					SYS(4, 0) = 0;
					SYS(4, 1) = y*(lam + 2 * mu);
					SYS(4, 2) = mu*x;
					SYS(4, 3) = 0;
					SYS(4, 4) = mu*x*x + y*y*(lam + 2 * mu);
					SYS(4, 5) = mu*x*z;
					SYS(4, 6) = y*z*(lam + 2 * mu);
					SYS(4, 7) = mu*x*x*z + y*y*z*(lam + 2 * mu);
					SYS(4, 8) = 0;
					SYS(4, 9) = mu*x;
					SYS(4, 10) = lam*y;
					SYS(4, 11) = 0;
					SYS(4, 12) = lam*x*y + mu*x*y;
					SYS(4, 13) = lam*y*z;
					SYS(4, 14) = mu*x*z;
					SYS(4, 15) = lam*x*y*z + mu*x*y*z;
					SYS(4, 16) = 0;
					SYS(4, 17) = 0;
					SYS(4, 18) = 0;
					SYS(4, 19) = lam*y;
					SYS(4, 20) = 0;
					SYS(4, 21) = lam*y*y;
					SYS(4, 22) = lam*x*y;
					SYS(4, 23) = lam*x*y*y;
					SYS(5, 0) = 0;
					SYS(5, 1) = 0;
					SYS(5, 2) = mu*z;
					SYS(5, 3) = mu*y;
					SYS(5, 4) = mu*x*z;
					SYS(5, 5) = mu*y*y + mu*z*z;
					SYS(5, 6) = mu*x*y;
					SYS(5, 7) = mu*x*y*y + mu*x*z*z;
					SYS(5, 8) = 0;
					SYS(5, 9) = mu*z;
					SYS(5, 10) = 0;
					SYS(5, 11) = 0;
					SYS(5, 12) = mu*y*z;
					SYS(5, 13) = 0;
					SYS(5, 14) = mu*z*z;
					SYS(5, 15) = mu*y*z*z;
					SYS(5, 16) = 0;
					SYS(5, 17) = mu*y;
					SYS(5, 18) = 0;
					SYS(5, 19) = 0;
					SYS(5, 20) = mu*y*y;
					SYS(5, 21) = 0;
					SYS(5, 22) = mu*y*z;
					SYS(5, 23) = mu*y*y*z;
					SYS(6, 0) = 0;
					SYS(6, 1) = z*(lam + 2 * mu);
					SYS(6, 2) = 0;
					SYS(6, 3) = mu*x;
					SYS(6, 4) = y*z*(lam + 2 * mu);
					SYS(6, 5) = mu*x*y;
					SYS(6, 6) = mu*x*x + z*z*(lam + 2 * mu);
					SYS(6, 7) = mu*x*x*y + y*z*z*(lam + 2 * mu);
					SYS(6, 8) = 0;
					SYS(6, 9) = 0;
					SYS(6, 10) = lam*z;
					SYS(6, 11) = 0;
					SYS(6, 12) = lam*x*z;
					SYS(6, 13) = lam*z*z;
					SYS(6, 14) = 0;
					SYS(6, 15) = lam*x*z*z;
					SYS(6, 16) = 0;
					SYS(6, 17) = mu*x;
					SYS(6, 18) = 0;
					SYS(6, 19) = lam*z;
					SYS(6, 20) = mu*x*y;
					SYS(6, 21) = lam*y*z;
					SYS(6, 22) = lam*x*z + mu*x*z;
					SYS(6, 23) = lam*x*y*z + mu*x*y*z;
					SYS(7, 0) = 0;
					SYS(7, 1) = y*z*(lam + 2 * mu);
					SYS(7, 2) = mu*x*z;
					SYS(7, 3) = mu*x*y;
					SYS(7, 4) = mu*x*x*z + y*y*z*(lam + 2 * mu);
					SYS(7, 5) = mu*x*y*y + mu*x*z*z;
					SYS(7, 6) = mu*x*x*y + y*z*z*(lam + 2 * mu);
					SYS(7, 7) = mu*x*x*y*y + mu*x*x*z*z + y*y*z*z*(lam + 2 * mu);
					SYS(7, 8) = 0;
					SYS(7, 9) = mu*x*z;
					SYS(7, 10) = lam*y*z;
					SYS(7, 11) = 0;
					SYS(7, 12) = lam*x*y*z + mu*x*y*z;
					SYS(7, 13) = lam*y*z*z;
					SYS(7, 14) = mu*x*z*z;
					SYS(7, 15) = lam*x*y*z*z + mu*x*y*z*z;
					SYS(7, 16) = 0;
					SYS(7, 17) = mu*x*y;
					SYS(7, 18) = 0;
					SYS(7, 19) = lam*y*z;
					SYS(7, 20) = mu*x*y*y;
					SYS(7, 21) = lam*y*y*z;
					SYS(7, 22) = lam*x*y*z + mu*x*y*z;
					SYS(7, 23) = lam*x*y*y*z + mu*x*y*y*z;
					SYS(8, 0) = 0;
					SYS(8, 1) = 0;
					SYS(8, 2) = 0;
					SYS(8, 3) = 0;
					SYS(8, 4) = 0;
					SYS(8, 5) = 0;
					SYS(8, 6) = 0;
					SYS(8, 7) = 0;
					SYS(8, 8) = 0;
					SYS(8, 9) = 0;
					SYS(8, 10) = 0;
					SYS(8, 11) = 0;
					SYS(8, 12) = 0;
					SYS(8, 13) = 0;
					SYS(8, 14) = 0;
					SYS(8, 15) = 0;
					SYS(8, 16) = 0;
					SYS(8, 17) = 0;
					SYS(8, 18) = 0;
					SYS(8, 19) = 0;
					SYS(8, 20) = 0;
					SYS(8, 21) = 0;
					SYS(8, 22) = 0;
					SYS(8, 23) = 0;
					SYS(9, 0) = 0;
					SYS(9, 1) = 0;
					SYS(9, 2) = mu;
					SYS(9, 3) = 0;
					SYS(9, 4) = mu*x;
					SYS(9, 5) = mu*z;
					SYS(9, 6) = 0;
					SYS(9, 7) = mu*x*z;
					SYS(9, 8) = 0;
					SYS(9, 9) = mu;
					SYS(9, 10) = 0;
					SYS(9, 11) = 0;
					SYS(9, 12) = mu*y;
					SYS(9, 13) = 0;
					SYS(9, 14) = mu*z;
					SYS(9, 15) = mu*y*z;
					SYS(9, 16) = 0;
					SYS(9, 17) = 0;
					SYS(9, 18) = 0;
					SYS(9, 19) = 0;
					SYS(9, 20) = 0;
					SYS(9, 21) = 0;
					SYS(9, 22) = 0;
					SYS(9, 23) = 0;
					SYS(10, 0) = 0;
					SYS(10, 1) = lam;
					SYS(10, 2) = 0;
					SYS(10, 3) = 0;
					SYS(10, 4) = lam*y;
					SYS(10, 5) = 0;
					SYS(10, 6) = lam*z;
					SYS(10, 7) = lam*y*z;
					SYS(10, 8) = 0;
					SYS(10, 9) = 0;
					SYS(10, 10) = lam + 2 * mu;
					SYS(10, 11) = 0;
					SYS(10, 12) = x*(lam + 2 * mu);
					SYS(10, 13) = z*(lam + 2 * mu);
					SYS(10, 14) = 0;
					SYS(10, 15) = x*z*(lam + 2 * mu);
					SYS(10, 16) = 0;
					SYS(10, 17) = 0;
					SYS(10, 18) = 0;
					SYS(10, 19) = lam;
					SYS(10, 20) = 0;
					SYS(10, 21) = lam*y;
					SYS(10, 22) = lam*x;
					SYS(10, 23) = lam*x*y;
					SYS(11, 0) = 0;
					SYS(11, 1) = 0;
					SYS(11, 2) = 0;
					SYS(11, 3) = 0;
					SYS(11, 4) = 0;
					SYS(11, 5) = 0;
					SYS(11, 6) = 0;
					SYS(11, 7) = 0;
					SYS(11, 8) = 0;
					SYS(11, 9) = 0;
					SYS(11, 10) = 0;
					SYS(11, 11) = mu;
					SYS(11, 12) = 0;
					SYS(11, 13) = mu*y;
					SYS(11, 14) = mu*x;
					SYS(11, 15) = mu*x*y;
					SYS(11, 16) = 0;
					SYS(11, 17) = 0;
					SYS(11, 18) = mu;
					SYS(11, 19) = 0;
					SYS(11, 20) = mu*x;
					SYS(11, 21) = mu*z;
					SYS(11, 22) = 0;
					SYS(11, 23) = mu*x*z;
					SYS(12, 0) = 0;
					SYS(12, 1) = lam*x;
					SYS(12, 2) = mu*y;
					SYS(12, 3) = 0;
					SYS(12, 4) = lam*x*y + mu*x*y;
					SYS(12, 5) = mu*y*z;
					SYS(12, 6) = lam*x*z;
					SYS(12, 7) = lam*x*y*z + mu*x*y*z;
					SYS(12, 8) = 0;
					SYS(12, 9) = mu*y;
					SYS(12, 10) = x*(lam + 2 * mu);
					SYS(12, 11) = 0;
					SYS(12, 12) = mu*y*y + x*x*(lam + 2 * mu);
					SYS(12, 13) = x*z*(lam + 2 * mu);
					SYS(12, 14) = mu*y*z;
					SYS(12, 15) = mu*y*y*z + x*x*z*(lam + 2 * mu);
					SYS(12, 16) = 0;
					SYS(12, 17) = 0;
					SYS(12, 18) = 0;
					SYS(12, 19) = lam*x;
					SYS(12, 20) = 0;
					SYS(12, 21) = lam*x*y;
					SYS(12, 22) = lam*x*x;
					SYS(12, 23) = lam*x*x*y;
					SYS(13, 0) = 0;
					SYS(13, 1) = lam*z;
					SYS(13, 2) = 0;
					SYS(13, 3) = 0;
					SYS(13, 4) = lam*y*z;
					SYS(13, 5) = 0;
					SYS(13, 6) = lam*z*z;
					SYS(13, 7) = lam*y*z*z;
					SYS(13, 8) = 0;
					SYS(13, 9) = 0;
					SYS(13, 10) = z*(lam + 2 * mu);
					SYS(13, 11) = mu*y;
					SYS(13, 12) = x*z*(lam + 2 * mu);
					SYS(13, 13) = mu*y*y + z*z*(lam + 2 * mu);
					SYS(13, 14) = mu*x*y;
					SYS(13, 15) = mu*x*y*y + x*z*z*(lam + 2 * mu);
					SYS(13, 16) = 0;
					SYS(13, 17) = 0;
					SYS(13, 18) = mu*y;
					SYS(13, 19) = lam*z;
					SYS(13, 20) = mu*x*y;
					SYS(13, 21) = lam*y*z + mu*y*z;
					SYS(13, 22) = lam*x*z;
					SYS(13, 23) = lam*x*y*z + mu*x*y*z;
					SYS(14, 0) = 0;
					SYS(14, 1) = 0;
					SYS(14, 2) = mu*z;
					SYS(14, 3) = 0;
					SYS(14, 4) = mu*x*z;
					SYS(14, 5) = mu*z*z;
					SYS(14, 6) = 0;
					SYS(14, 7) = mu*x*z*z;
					SYS(14, 8) = 0;
					SYS(14, 9) = mu*z;
					SYS(14, 10) = 0;
					SYS(14, 11) = mu*x;
					SYS(14, 12) = mu*y*z;
					SYS(14, 13) = mu*x*y;
					SYS(14, 14) = mu*x*x + mu*z*z;
					SYS(14, 15) = mu*x*x*y + mu*y*z*z;
					SYS(14, 16) = 0;
					SYS(14, 17) = 0;
					SYS(14, 18) = mu*x;
					SYS(14, 19) = 0;
					SYS(14, 20) = mu*x*x;
					SYS(14, 21) = mu*x*z;
					SYS(14, 22) = 0;
					SYS(14, 23) = mu*x*x*z;
					SYS(15, 0) = 0;
					SYS(15, 1) = lam*x*z;
					SYS(15, 2) = mu*y*z;
					SYS(15, 3) = 0;
					SYS(15, 4) = lam*x*y*z + mu*x*y*z;
					SYS(15, 5) = mu*y*z*z;
					SYS(15, 6) = lam*x*z*z;
					SYS(15, 7) = lam*x*y*z*z + mu*x*y*z*z;
					SYS(15, 8) = 0;
					SYS(15, 9) = mu*y*z;
					SYS(15, 10) = x*z*(lam + 2 * mu);
					SYS(15, 11) = mu*x*y;
					SYS(15, 12) = mu*y*y*z + x*x*z*(lam + 2 * mu);
					SYS(15, 13) = mu*x*y*y + x*z*z*(lam + 2 * mu);
					SYS(15, 14) = mu*x*x*y + mu*y*z*z;
					SYS(15, 15) = mu*x*x*y*y + mu*y*y*z*z + x*x*z*z*(lam + 2 * mu);
					SYS(15, 16) = 0;
					SYS(15, 17) = 0;
					SYS(15, 18) = mu*x*y;
					SYS(15, 19) = lam*x*z;
					SYS(15, 20) = mu*x*x*y;
					SYS(15, 21) = lam*x*y*z + mu*x*y*z;
					SYS(15, 22) = lam*x*x*z;
					SYS(15, 23) = lam*x*x*y*z + mu*x*x*y*z;
					SYS(16, 0) = 0;
					SYS(16, 1) = 0;
					SYS(16, 2) = 0;
					SYS(16, 3) = 0;
					SYS(16, 4) = 0;
					SYS(16, 5) = 0;
					SYS(16, 6) = 0;
					SYS(16, 7) = 0;
					SYS(16, 8) = 0;
					SYS(16, 9) = 0;
					SYS(16, 10) = 0;
					SYS(16, 11) = 0;
					SYS(16, 12) = 0;
					SYS(16, 13) = 0;
					SYS(16, 14) = 0;
					SYS(16, 15) = 0;
					SYS(16, 16) = 0;
					SYS(16, 17) = 0;
					SYS(16, 18) = 0;
					SYS(16, 19) = 0;
					SYS(16, 20) = 0;
					SYS(16, 21) = 0;
					SYS(16, 22) = 0;
					SYS(16, 23) = 0;
					SYS(17, 0) = 0;
					SYS(17, 1) = 0;
					SYS(17, 2) = 0;
					SYS(17, 3) = mu;
					SYS(17, 4) = 0;
					SYS(17, 5) = mu*y;
					SYS(17, 6) = mu*x;
					SYS(17, 7) = mu*x*y;
					SYS(17, 8) = 0;
					SYS(17, 9) = 0;
					SYS(17, 10) = 0;
					SYS(17, 11) = 0;
					SYS(17, 12) = 0;
					SYS(17, 13) = 0;
					SYS(17, 14) = 0;
					SYS(17, 15) = 0;
					SYS(17, 16) = 0;
					SYS(17, 17) = mu;
					SYS(17, 18) = 0;
					SYS(17, 19) = 0;
					SYS(17, 20) = mu*y;
					SYS(17, 21) = 0;
					SYS(17, 22) = mu*z;
					SYS(17, 23) = mu*y*z;
					SYS(18, 0) = 0;
					SYS(18, 1) = 0;
					SYS(18, 2) = 0;
					SYS(18, 3) = 0;
					SYS(18, 4) = 0;
					SYS(18, 5) = 0;
					SYS(18, 6) = 0;
					SYS(18, 7) = 0;
					SYS(18, 8) = 0;
					SYS(18, 9) = 0;
					SYS(18, 10) = 0;
					SYS(18, 11) = mu;
					SYS(18, 12) = 0;
					SYS(18, 13) = mu*y;
					SYS(18, 14) = mu*x;
					SYS(18, 15) = mu*x*y;
					SYS(18, 16) = 0;
					SYS(18, 17) = 0;
					SYS(18, 18) = mu;
					SYS(18, 19) = 0;
					SYS(18, 20) = mu*x;
					SYS(18, 21) = mu*z;
					SYS(18, 22) = 0;
					SYS(18, 23) = mu*x*z;
					SYS(19, 0) = 0;
					SYS(19, 1) = lam;
					SYS(19, 2) = 0;
					SYS(19, 3) = 0;
					SYS(19, 4) = lam*y;
					SYS(19, 5) = 0;
					SYS(19, 6) = lam*z;
					SYS(19, 7) = lam*y*z;
					SYS(19, 8) = 0;
					SYS(19, 9) = 0;
					SYS(19, 10) = lam;
					SYS(19, 11) = 0;
					SYS(19, 12) = lam*x;
					SYS(19, 13) = lam*z;
					SYS(19, 14) = 0;
					SYS(19, 15) = lam*x*z;
					SYS(19, 16) = 0;
					SYS(19, 17) = 0;
					SYS(19, 18) = 0;
					SYS(19, 19) = lam + 2 * mu;
					SYS(19, 20) = 0;
					SYS(19, 21) = y*(lam + 2 * mu);
					SYS(19, 22) = x*(lam + 2 * mu);
					SYS(19, 23) = x*y*(lam + 2 * mu);
					SYS(20, 0) = 0;
					SYS(20, 1) = 0;
					SYS(20, 2) = 0;
					SYS(20, 3) = mu*y;
					SYS(20, 4) = 0;
					SYS(20, 5) = mu*y*y;
					SYS(20, 6) = mu*x*y;
					SYS(20, 7) = mu*x*y*y;
					SYS(20, 8) = 0;
					SYS(20, 9) = 0;
					SYS(20, 10) = 0;
					SYS(20, 11) = mu*x;
					SYS(20, 12) = 0;
					SYS(20, 13) = mu*x*y;
					SYS(20, 14) = mu*x*x;
					SYS(20, 15) = mu*x*x*y;
					SYS(20, 16) = 0;
					SYS(20, 17) = mu*y;
					SYS(20, 18) = mu*x;
					SYS(20, 19) = 0;
					SYS(20, 20) = mu*x*x + mu*y*y;
					SYS(20, 21) = mu*x*z;
					SYS(20, 22) = mu*y*z;
					SYS(20, 23) = mu*x*x*z + mu*y*y*z;
					SYS(21, 0) = 0;
					SYS(21, 1) = lam*y;
					SYS(21, 2) = 0;
					SYS(21, 3) = 0;
					SYS(21, 4) = lam*y*y;
					SYS(21, 5) = 0;
					SYS(21, 6) = lam*y*z;
					SYS(21, 7) = lam*y*y*z;
					SYS(21, 8) = 0;
					SYS(21, 9) = 0;
					SYS(21, 10) = lam*y;
					SYS(21, 11) = mu*z;
					SYS(21, 12) = lam*x*y;
					SYS(21, 13) = lam*y*z + mu*y*z;
					SYS(21, 14) = mu*x*z;
					SYS(21, 15) = lam*x*y*z + mu*x*y*z;
					SYS(21, 16) = 0;
					SYS(21, 17) = 0;
					SYS(21, 18) = mu*z;
					SYS(21, 19) = y*(lam + 2 * mu);
					SYS(21, 20) = mu*x*z;
					SYS(21, 21) = mu*z*z + y*y*(lam + 2 * mu);
					SYS(21, 22) = x*y*(lam + 2 * mu);
					SYS(21, 23) = mu*x*z*z + x*y*y*(lam + 2 * mu);
					SYS(22, 0) = 0;
					SYS(22, 1) = lam*x;
					SYS(22, 2) = 0;
					SYS(22, 3) = mu*z;
					SYS(22, 4) = lam*x*y;
					SYS(22, 5) = mu*y*z;
					SYS(22, 6) = lam*x*z + mu*x*z;
					SYS(22, 7) = lam*x*y*z + mu*x*y*z;
					SYS(22, 8) = 0;
					SYS(22, 9) = 0;
					SYS(22, 10) = lam*x;
					SYS(22, 11) = 0;
					SYS(22, 12) = lam*x*x;
					SYS(22, 13) = lam*x*z;
					SYS(22, 14) = 0;
					SYS(22, 15) = lam*x*x*z;
					SYS(22, 16) = 0;
					SYS(22, 17) = mu*z;
					SYS(22, 18) = 0;
					SYS(22, 19) = x*(lam + 2 * mu);
					SYS(22, 20) = mu*y*z;
					SYS(22, 21) = x*y*(lam + 2 * mu);
					SYS(22, 22) = mu*z*z + x*x*(lam + 2 * mu);
					SYS(22, 23) = mu*y*z*z + x*x*y*(lam + 2 * mu);
					SYS(23, 0) = 0;
					SYS(23, 1) = lam*x*y;
					SYS(23, 2) = 0;
					SYS(23, 3) = mu*y*z;
					SYS(23, 4) = lam*x*y*y;
					SYS(23, 5) = mu*y*y*z;
					SYS(23, 6) = lam*x*y*z + mu*x*y*z;
					SYS(23, 7) = lam*x*y*y*z + mu*x*y*y*z;
					SYS(23, 8) = 0;
					SYS(23, 9) = 0;
					SYS(23, 10) = lam*x*y;
					SYS(23, 11) = mu*x*z;
					SYS(23, 12) = lam*x*x*y;
					SYS(23, 13) = lam*x*y*z + mu*x*y*z;
					SYS(23, 14) = mu*x*x*z;
					SYS(23, 15) = lam*x*x*y*z + mu*x*x*y*z;
					SYS(23, 16) = 0;
					SYS(23, 17) = mu*y*z;
					SYS(23, 18) = mu*x*z;
					SYS(23, 19) = x*y*(lam + 2 * mu);
					SYS(23, 20) = mu*x*x*z + mu*y*y*z;
					SYS(23, 21) = mu*x*z*z + x*y*y*(lam + 2 * mu);
					SYS(23, 22) = mu*y*z*z + x*x*y*(lam + 2 * mu);
					SYS(23, 23) = mu*x*x*z*z + mu*y*y*z*z + x*x*y*y*(lam + 2 * mu);
					//Ïåðåîïðåäåëåíèå SYS>

					for (u1 = 0; u1 < q; u1++)
						for (u2 = 0; u2 < q; u2++)
							K(u1, u2) = K(u1, u2) + SYS(u1, u2) * CTData[j + m2*k + s*h][i];

				}

		double VOX = voxX*voxY*voxZ;

		K = VOX * K;
		//Èíòåãðèðîâàíèå ïî áèíàðèçîâàííîìó ìàññèâó>

		//<Ïîñòðîåíèå À è AT
		for (i = 0; i < q; i++)
			for (j = 0; j < q; j++)
				A(i, j) = 0;


		/*for (j = 0; j < q; j += n) {

		for (i = 0; i < n; i++)

		A(i + j, j) = 1;

		for (i = 1; i <= 2; i++) {
		for (int i0 = 0; i0 <= 4; i0 += 4) {

		A(i + j + i0, j + 1) = a;
		A(i + j + i0 + 1, j + 2) = b;
		A(i + j + i0 / 2 + 3, j + 3) = c;

		A(i + j - i0 / 4 + 5, j + i0 / 4 + 5) = b*c;

		}

		A((i - 1) * 4 + j + 2, j + 4) = a*b;

		}

		A(6 + j, j + 7) = a*b*c;

		}*/


		/*Êàê îïðåäåëèòü êîîðäíàòû óçëà? Ëåãêî! Äîñòàòî÷íî ïðèäåðæèâàòüñÿ ñëåäóþùèõ èíñòðóêöèé:

		Nodes[ Elems[a1][a2] ][ a3 ]:
		a1 - íîìåð ýëåìåíòà,
		a2 - íîìåð óçëà (ëîêàëüíûé),
		a3 - íîìåð êîîðäèíàòû (x, y èëè z).

		Ïóñü a=(a1,a2,a3), òîãäà a=(0,0,0) ñîîòâåòñòâóåò x-êîîðäèíàòà 0-ãî óçëà 0-ãî ýëåìåíòà.*/

		Eigen::VectorXd CLx(n), CLy(n), CLz(n);	//Ëîêàëüíûå êîîðäèíàòû (äëÿ îïðåäåëåíèÿ ìàòðèöû A).

												//<Ðàçìåð òåêóùåãî ÊÝ
		a = Nodes[Elems[s][n - 1]][0] - Nodes[Elems[s][0]][0];
		b = Nodes[Elems[s][n - 1]][1] - Nodes[Elems[s][0]][1];
		c = Nodes[Elems[s][n - 1]][2] - Nodes[Elems[s][0]][2];
		//Ðàçìåð òåêóùåãî ÊÝ>

		//<Ëîêàëüíûå êîîðäèíàòû òåêóùåãî ÊÝ
		for (i = 0; i < 2; i++) {

			for (j = 0; j < 4; j++)
				CLx(i + 2 * j) = i*a;

			for (j = 0; j <= 4; j += 4)
				for (k = 0; k < 2; k++)
					CLy(2 * i + j + k) = i*b;

			for (j = 0; j < 4; j++)
				CLz(4 * i + j) = i*c;

		}
		//Ëîêàëüíûå êîîðäèíàòû òåêóùåãî ÊÝ>

		for (i = 0; i < n; i++) {

			cout << CLx(i) << ' ' << CLy(i) << ' ' << CLz(i) << endl;

			for (j = 0; j <= 2 * n; j += n) {

				A(i + j, j) = 1;
				A(i + j, j + 1) = CLx(i);
				A(i + j, j + 2) = CLy(i);
				A(i + j, j + 3) = CLz(i);
				A(i + j, j + 4) = CLx(i)*CLy(i);
				A(i + j, j + 5) = CLy(i)*CLz(i);
				A(i + j, j + 6) = CLx(i)*CLz(i);
				A(i + j, j + 7) = CLx(i)*CLy(i)*CLz(i);

			}
		}

		AT = A.inverse();
		//Ïîñòðîåíèå À è AT>

		//	//<À è ÀÒ
		std::ofstream outfileA("C:/Users/imm/Documents/Visual Studio 2017/Projects/ConsoleApplication1/A.txt", std::ofstream::out);
		outfileA << A << std::endl;
		outfileA << std::endl;
		outfileA << AT << std::endl;
		//	outfileA.close();

		//<Ïîñòðîåíèå ËÌÆ
		KL = AT.transpose()*K*AT;
		//Ïîñòðîåíèå ËÌÆ>

	}

public:
	//Ïîñòðîåíèå ÃÌÆ:
	void GlobalMatrix() {

		SparseKG.resize(NNodes*p, NNodes*p);

		int c1, c2;
		int d1, d2, d3, d4;

		for (s = 0; s < NElems; s++) {

			/*Êàê ëåãêî îïðåäåëèòü ãëîáàëüíûé íîìåð óçëà? Äîñòàòî÷íî èñïîëüçîâàòü ñëåäóþùóþ èíñòðóêöèþ:

			Elems[b1][b2]:
			b1 - íîìåð ýëåìåíòà,
			b2 - íîìåð óçëà (ëîêàëüíûé).*/

			this->LocalMatrix();

			//<Ñóììèðîâàíèå
			for (i = 0; i < n; i++)
				for (j = 0; j < n; j++) {

					c1 = Elems[s][i];
					c2 = Elems[s][j];

					for (int k1 = 0; k1 <= 2; k1++)
						for (int k2 = 0; k2 <= 2; k2++) {

							d1 = c1 + k1*NNodes;
							d2 = c2 + k2*NNodes;

							d3 = i + k1*n;
							d4 = j + k2*n;

							SparseKG.coeffRef(d1, d2) = SparseKG.coeff(d1, d2) + KL(d3, d4);

						}
				}
			//Ñóììèðîâàíèå>

		}

		//<KG
		std::ofstream outfileSparseKG("C:/Users/imm/Documents/Visual Studio 2017/Projects/ConsoleApplication1/SparseKG.txt", std::ofstream::out);
		outfileSparseKG << SparseKG << std::endl;
		outfileSparseKG << std::endl;
		outfileSparseKG.close();
		//KG>

	}

	//Ðåøåíèå â ïåðåìåùåíèÿõ:
	void SolveDelta() {

		this->GlobalMatrix();

		//<"Èñêëþ÷åíèå" ñòðîê è ñòîëáöîâ
		for (i = 0; i < SparseKG.outerSize(); i++)
			for (Eigen::SparseMatrix<double>::InnerIterator it(SparseKG, i); it; ++it)
				for (std::vector<int>::iterator idit = BCKinem.begin(); idit != BCKinem.end(); ++idit)
					if (it.row() == *idit || it.col() == *idit)
						//if (it.row() == *idit && it.col() == *idit) it.valueRef() = 1.0/**it.value()*/; else it.valueRef() = 0.0;
						it.valueRef() = it.row() == it.col() ? 1.0 : 0.0;
		//"Èñêëþ÷åíèå" ñòðîê è ñòîëáöîâ>

		//<Ðåøåíèå ÑËÀÓ
		Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > solver(SparseKG);
		Delta = solver.solve(BCForce);
		//Ðåøåíèå ÑËÀÓ>

		std::cout << Delta << std::endl;

	}

	//Çàïèñü â ôàéë:
	//void WriteToFile(Eigen::MatrixXd &K, Eigen::MatrixXd &A, Eigen::MatrixXd &AT, 
	//	Eigen::MatrixXd &KL, Eigen::SparseMatrix<double> &SparseKL, Eigen::VectorXd &Delta) {
	//
	//	//<K
	//	std::ofstream outfileK("C:/Users/imm/Documents/Visual Studio 2010/Projects/1/K.txt", std::ofstream::out);
	//	outfileK << K << std::endl;
	//	outfileK.close();
	//	//K>
	//
	//	//<À è ÀÒ
	//	std::ofstream outfileA("C:/Users/imm/Documents/Visual Studio 2010/Projects/1/A.txt", std::ofstream::out);
	//	outfileA << A << std::endl;
	//	outfileA << std::endl;
	//	outfileA << AT << std::endl;
	//	outfileA.close();
	//
	//	//Âûâîä íà ìîíèòîð:
	//	//std::cout << A << endl;
	//	//std::cout << AT << endl;
	//	//À è ÀÒ>
	//
	//	//<KL
	//	std::ofstream outfileKL("C:/Users/imm/Documents/Visual Studio 2010/Projects/1/KL.txt", std::ofstream::out);
	//	outfileKL << KL << std::endl;
	//	outfileKL.close();
	//	//KL>
	//
	//	//<SparseKL
	//	std::ofstream outfileSparseKL("C:/Users/imm/Documents/Visual Studio 2010/Projects/1/SparseKL.txt", std::ofstream::out);
	//	outfileSparseKL << SparseKL << std::endl;
	//	outfileSparseKL.close();
	//	//SparseKL>
	//
	//	//<Delta
	//	std::ofstream outfileDelta("C:/Users/imm/Documents/Visual Studio 2010/Projects/1/Delta.txt", std::ofstream::out);
	//	outfileDelta << Delta << std::endl;
	//	outfileDelta.close();
	//	//Delta>
	//
	//}

};

int main() {

	setlocale(LC_ALL, "RUSSIAN");

	gFEM3D8Node* Solve = new gFEM3D8Node();

	Solve->SolveDelta();
	//Solve->GlobalMatrix();

	return 0;
}

