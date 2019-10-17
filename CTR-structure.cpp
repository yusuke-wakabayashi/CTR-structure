#include<fstream>
#include<iostream>
#include<vector>
#include<string>
#include<math.h>
#include<random>
#include<complex>
#include<stdio.h>
#include<array>
#include<time.h>
#include<direct.h>
#define I (std::complex<double> (0,1))
#define PI 3.1415926535
// Scale of the random walk for parameter search パラメーター探索ランダムウォークの移動幅
#define DELTA 0.01
// MCMC sampling number between the temperature exchange  １交換間に行うMCMCサンプリング数
#define LOOP_NUM 100
//  MCMC sampling number between the output of .par file.  .parファイル出力頻度
#define RECORD_NUM 10
//  Virtual linear absorption coefficient (in the unit of /unit cell).  X線吸収率(/u.c.)
#define MU 0.02 
// Number of replicas for exchange Monte Carlo calc.   交換モンテカルロ法の平行計算レプリカ数
#define REPLICA_NUM 128


// Ignore comments written in Japanese. They are not very important, but left for future maintenance.

namespace CTRstructure {//グローバル空間汚染防止
	//標準C++ライブラリのうち以下はstdを用いずに使用する
	using std::array;
	using std::vector;

	using std::complex;
	using std::ifstream;
	using std::ofstream;
	using std::string;

	using std::exception;
	using std::runtime_error;
	using std::invalid_argument;

	using std::to_string;

	//使用するクラス一覧（ヘッダー替わり）
	class F0Data;
	class AtomData;
	class BulkData;
	class FilmData;
	class ExpData;
	class Replica;

	//共通で使用するプロパティ
	//ログ出力用のファイル
	ofstream logFile("log.txt", std::ios::app);

	//乱数は毎回異なる初期化シードを用いる
	// Random number seed is changed for each run.
	std::random_device seed_gen;
	std::mt19937_64 random(seed_gen());

	/*汎用関数　現在時刻を取得する*/
	// General function: obtain the current time.
	string getTime() {
		char date[64];
		time_t t = time(NULL);
		struct tm tm;
		localtime_s(&tm, &t);
		strftime(date, sizeof(date), "%Y/%m/%d %a %H:%M:%S", &tm);
		string result(date);
		return result;
	}

	/*
	汎用関数
	ログファイルにstringを書き出す
	General function: output string to the log file.
	*/
	void outLog(const string str) {
		logFile << str << std::endl;
		std::cout << str << std::endl;
	}
	void outError(const exception& e, const int line) {
		logFile << "ERROR:" << e.what() << ",line:" << line << std::endl;
		std::cout << "ERROR:" << e.what() << ",line:" << line << std::endl;
	}
	/*
	汎用関数
	文字列strを文字regexにて分割し、vector形式の配列で返す
	javaのstr.split(String regex,-1)と同じ処理
	strに2バイト文字を入れた場合どうなるかは不明
	Gneneral function: split a string.
	*/
	vector<string> splitStr(const string str, const char regex) {
		try {
			vector<string> elems;
			string item;
			for (char ch : str) {
				if (ch == regex) {
					elems.push_back(item);
					item.clear();
				}
				else {
					item += ch;
				}
			}
			if (!item.empty())
				elems.push_back(item);
			return elems;
		}
		catch (exception& e) {
			outError(e, __LINE__);
			throw e;
		}
		catch (...) {
			runtime_error e("Unknown error");
			outError(e, __LINE__);
			throw runtime_error("Unknown error");
		}
	}
	/**
	乱数生成器
	gaussSwitchをTrueにするとガウス乱数、falseにすると０から１の一様乱数
	Random number generator.
	gaussSwitch=True  :  Gaussian random numbers,
	gaussSwitch=false :  uniform random numbers,
	*/
	double gauss_rand(const bool gaussSwitch) {
		double r1 = 1.0 - (double)random() / ((double)random._Max + 1);// divergence stopper for log calc. log内の発散防止
		double r2 = 1.0 - (double)random() / ((double)random._Max + 1);
		double r3;

		if (gaussSwitch)
			return sqrt(-2.0 * log(r1)) * sin(2.0 * PI * r2);// Gauss rand with average of 0 and standard deviation of 1.  平均0,分散1の正規乱数。パラメータを動かす用

		else {
			//パラメータ判定用。
			r3 = (double)random() / ((double)random._Max);
			return  r3;// uniform random number between 0 and 1.  0から1の一様乱数
		}
	}

	/*
	原子散乱因子のデータを管理する
	異常分散構造のデータもここで管理する
	Atomic scattering factor and anomalous dispersion are treated in this class.
	*/
	class F0Data {
	public:
		static const int MAX_ELEMENTS = 150;

		//これらの配列はH=0,He=1の順で格納されていることに注意
		// The first argument corresponds to the value (atomic number -1). e.g., fp[1] has the value of f' for helium, not for hydrogen.
		double A[MAX_ELEMENTS][4] = { 0 };
		double B[MAX_ELEMENTS][4] = { 0 };
		double C[MAX_ELEMENTS] = { 0 };
		double fp[MAX_ELEMENTS] = { 0 };
		double fpp[MAX_ELEMENTS] = { 0 };
		F0Data(const string inputFileName) {
			try {
				ifstream inputFile(inputFileName);
				if (!inputFile) {
					throw runtime_error(string("FileNotFoundException(") + inputFileName + string(")"));
				}
				//読み込み時のみに使う変数
				// temporal parameters for read out.
				string str;
				int atomicNum = 0;
				while (getline(inputFile, str)) {
					if (atomicNum < MAX_ELEMENTS) {
						vector<string> strList = splitStr(str, ' ');
						if (strList.size() >= 9) {
							for (size_t i = 0; i < 4; i++) {
								A[atomicNum][i] = stod(strList[2 * i]);
								B[atomicNum][i] = stod(strList[2 * i + 1]);
							}
							C[atomicNum] = stod(strList[8]);
							atomicNum++;
						}
					}
					else {
						throw runtime_error(string("Too long f0data.dat. MAX_ELEMENTS=" + to_string(MAX_ELEMENTS) + string(".")));
					}
				}
			}
			catch (exception& e) {
				outError(e, __LINE__);
				throw e;
			}
			catch (...) {
				runtime_error e("Unknown error");
				outError(e, __LINE__);
				throw e;
			}
		}
		int anomalousLoad(const string fileName) {
			try {
				ifstream inputFile(fileName);
				if (!inputFile) {
					outLog("No anomalousFile was found. The values of fp and fpp are set to 0.");
					return -1;
				}

				string str;
				while (getline(inputFile, str)) {
					vector<string> strList = splitStr(str, '\t');
					if (strList.size() >= 3) {
						fp[stoi(strList[0])] = stod(strList[1]);
						fpp[stoi(strList[0])] = stod(strList[2]);
					}
				}
				return 0;
			}
			catch (exception& e) {
				outError(e, __LINE__);
				throw e;
			}
			catch (...) {
				runtime_error e("Unknown error");
				outError(e, __LINE__);
				throw e;
			}


		}

		complex<double> f(const int atomicNumber, const double q)  const {
			try {
				if (atomicNumber <= 0) {
					return 0;
				}
				if (atomicNumber > MAX_ELEMENTS) {
					runtime_error e("ArrayIndexOutOfBoundsException(atomicNumber >MAX_ELEMENTS)");
					outError(e, __LINE__);
					throw e;
				}
				double result = C[atomicNumber - 1];
				for (int i = 0; i < 4; i++) {
					result += A[atomicNumber - 1][i] * exp(-B[atomicNumber - 1][i] * q * q);
				}
				return result + fp[atomicNumber] + I * fpp[atomicNumber];
			}
			catch (exception& e) {
				outError(e, __LINE__);
				throw e;
			}
			catch (...) {
				runtime_error e("Unknown error");
				outError(e, __LINE__);
				throw e;
			}
		}
	};


	class AtomData {
	public:

		int atomicNumber = 0;
		int atomSum = 0;
		static const int PARAM_NUM = 7;
		/*
		 index of param
		*/
		static const int d = 0;
		static const int x = 1;
		static const int y = 2;
		static const int z = 3;
		static const int B = 4;
		static const int toccf = 5;
		static const int foccf = 6;

		array<double, PARAM_NUM>param = { 0 };
		array<bool, PARAM_NUM> move = { false };
		array<int, PARAM_NUM>connectNum = { 0 };
		array<vector<int>, PARAM_NUM>connectAtomIndex;
		array<vector<int>, PARAM_NUM>connectParIndex;
		array<vector<double>, PARAM_NUM>connectRatio;

		array<double, PARAM_NUM> prospectMin = { 0 };
		array<double, PARAM_NUM> prospectMax = { 0 };
		complex<double> calcF(const double h, const double k, const double l, const double q, const complex<double>f) const {
			if (atomSum == 1) {
				return f * param[toccf] * param[foccf]
					* exp(2.0 * PI * I * (h * param[x] + k * param[y] + l * param[z]))
					* exp(-param[B] * q * q) * exp(MU * param[z]);
			}
			else {
				complex<double> result = 0;
				for (int i = 0; i < atomSum; i++) {
					result += f * param[toccf] * param[foccf]
						* exp(2.0 * PI * I * (h * param[x] + k * param[y] + l * (param[z] + i * param[d])))
						* exp(-param[B] * q * q) * exp(MU * (param[z] + i * param[d]));
				}
				return result;
			}
		}
		void setAtomicNumber(const int value) {
			if (value >= F0Data::MAX_ELEMENTS) {
				throw runtime_error("Too large atomicNumber");
			}
			atomicNumber = value;
		}

		static string getParamName(const int value) {
			switch (value)
			{
			case 0:
				return("d");
			case 1:
				return("x");
			case 2:
				return("y");
			case 3:
				return("z");
			case 4:
				return("B");
			case 5:
				return("toccf");
			case 6:
				return("foccf");
			}
			return "";
		}

		void addConnect(const int connectParam, const int atomIndex, const int paramIndex, const string ratio) {
			try {
				connectNum.at(connectParam)++;
				connectAtomIndex.at(connectParam).push_back(atomIndex);
				connectParIndex.at(connectParam).push_back(paramIndex);
				connectRatio.at(connectParam).push_back(stod(ratio));
			}
			catch (std::invalid_argument& e) {
				outLog("can't cast '" + ratio + "' to double value");
				outError(e, __LINE__);
				throw e;
			}
			catch (std::exception& e) {
				outError(e, __LINE__);
				throw e;
			}
		}




	};
	class FilmData {
	public:
		static const int MAX_ATOM_NUM = 100;

		double scale = 0;
		double noiseScale = 0;
		double noiseValue = 0;
		double bulkB = 0;
		bool scaleFit = false;
		bool noiseScaleFit = false;
		bool noiseValueFit = false;
		bool bulkBfit = false;

		array<AtomData, MAX_ATOM_NUM> atomData;
		int atomNum = 0;

		FilmData() {}
		FilmData(const string inputFileName) {
			ifstream inputFile(inputFileName);
			if (!inputFile) {
				throw runtime_error(string("FileNotFoundException(") + inputFileName + string(")"));
			}
			string str1;//読み込み用変数
			try {
				getline(inputFile, str1);
				vector<string> strList1 = splitStr(str1, '\t');
				scale = stod(strList1[1]);
				scaleFit = stoi(strList1[3]);
				noiseScale = stod(strList1[5]);
				noiseScaleFit = stoi(strList1[7]);
				noiseValue = stod(strList1[9]);
				noiseValueFit = stoi(strList1[11]);
				if (strList1.size() > 15) {
					bulkB = stod(strList1[13]);
					bulkBfit = stoi(strList1[15]);
				}
				else {
					bulkB = -1;
					bulkBfit = false;
				}

				getline(inputFile, str1);//捨て行
			}
			catch (...) {
				throw runtime_error(string("first line of the  film file->") + str1);
			}
			string str;//読み込み変数
			while (getline(inputFile, str)) {
				try {
					if (atomNum >= MAX_ATOM_NUM) {
						throw runtime_error(string("Too many atoms in the film file ->" + to_string(atomNum + 1) + string("MAX_ATOM_NUM")));
					}
					vector<string>strList = splitStr(str, '\t');
					strList.resize(50);
					atomData[atomNum].setAtomicNumber(stoi(strList[1]));
					atomData[atomNum].atomSum = stoi(strList[2]);
					for (int i = 0; i < AtomData::PARAM_NUM; i++) {
						atomData[atomNum].param[i] = stod(strList[i + 3]);
					}
					int fitPar = stoi(strList[10]);
					for (int i = 0; i < AtomData::PARAM_NUM; i++) {
						atomData[atomNum].move[AtomData::PARAM_NUM - 1 - i] = fitPar % 10;
						fitPar = fitPar / 10;
					}
					//接続パラメーター
					for (int i = 0; i < AtomData::PARAM_NUM; i++) {
						vector<string>connectAtomStr = splitStr(strList[3 * i + 11], ',');
						vector<string>connectParStr = splitStr(strList[3 * i + 12], ',');
						vector<string>connectRatioStr = splitStr(strList[3 * i + 13], ',');
						if ((connectAtomStr.size() != connectParStr.size()) || (connectAtomStr.size() != connectRatioStr.size())) {
							throw runtime_error(string("The number of connection conditions of atoms must be the same. Check ") + AtomData::getParamName(i) + string(" paramter's connection"));
						}
						for (int j = 0; j < connectAtomStr.size(); j++) {
							try {
								int connectAtom = stoi(connectAtomStr[j]);
								int connectPar = stoi(connectParStr[j]);
								atomData.at(connectAtom).addConnect(connectPar, atomNum, i, connectRatioStr[j]);
							}
							catch (exception& e) {
								outLog(string("Wrong connection info. Check ") + AtomData::getParamName(i) + string(" paramter's connection"));
								throw e;
							}
						}
					}
					//拘束パラメーター
					for (long i = 0; i < AtomData::PARAM_NUM; i++) {
						try {
							atomData[atomNum].prospectMin[i] = stod(strList[2 * i + 32]);
						}
						catch (invalid_argument e) {
							outLog(string("cast error in prospectMin. Software use the default value -10000"));
							atomData[atomNum].prospectMin[i] = -10000;
						}
						try {
							atomData[atomNum].prospectMax[i] = stod(strList[2 * i + 33]);
						}
						catch (invalid_argument e) {
							outLog(string("cast error in prospectMax. Software use the default value 10000"));
							atomData[atomNum].prospectMax[i] = 10000;
						}
					}
				}
				catch (exception& e) {
					outLog(string("Check the line:") + str);
					outError(e, __LINE__);
					throw e;
				}
				catch (...) {
					outLog(string("Check the line:") + str);
					runtime_error e("Unknown error");
					outError(e, __LINE__);
					throw e;
				}
				atomNum++;
			}
			inputFile.close();
		}
		const void outPosFile(const string outName)const {
			array<array<vector<int>, AtomData::PARAM_NUM>, MAX_ATOM_NUM>tmpIndex;
			array<array<vector<int>, AtomData::PARAM_NUM>, MAX_ATOM_NUM>tmpParam;
			array<array<vector<double>, AtomData::PARAM_NUM>, MAX_ATOM_NUM>tmpRatio;
			array<array<int, AtomData::PARAM_NUM>, MAX_ATOM_NUM>tmpNum = { 0 };

			try {

				for (int n = 0; n < atomNum; n++) {
					for (int i = 0; i < AtomData::PARAM_NUM; i++) {
						for (int k = 0; k < atomData[n].connectNum[i]; k++) {
							tmpNum.at(atomData[n].connectAtomIndex[i][k]).at(atomData[n].connectParIndex[i][k])++;
							tmpIndex.at(atomData[n].connectAtomIndex[i][k]).at(atomData[n].connectParIndex[i][k]).push_back(n);
							tmpParam.at(atomData[n].connectAtomIndex[i][k]).at(atomData[n].connectParIndex[i][k]).push_back(i);
							tmpRatio.at(atomData[n].connectAtomIndex[i][k]).at(atomData[n].connectParIndex[i][k]).push_back(atomData[n].connectRatio[i][k]);
						}
					}
				}

				ofstream outFile(outName);
				outFile << "scale\t" << scale << "\tscale_switch\t" << scaleFit << "\tgamma\t" << noiseScale << "\tgamma_switch\t" << noiseScaleFit << "\tsigma\t" << noiseValue << "\tsigma_switch\t" << noiseValueFit << "\tsubstrate-B\t" << bulkB << "\tsubstrate-B_switch\t" << bulkBfit << std::endl;
				outFile << "AtomIndex	AtmNum	N	d	x	y	z	B	toccf	foccf	7dig_fit_switch	d-index	d-param	d-ratio	x-index	x-param	x-ratio	y-index	y-param	y-ratio	z-index	z-param	z-ratio	B-index	B-param	B-ratio	toccf-index	toccf-param	toccf-ratio	foccf-index	foccf-param	foccf-ratio	d-min	d-max	x-min	x-max	y-min	y-max	z-min	z-max	B-min	B-max	toccf-min	toccf-max	foccf-min	foccf-max" << std::endl;
				for (int n = 0; n < atomNum; n++) {
					outFile << n << "\t";
					outFile << atomData[n].atomicNumber << "\t";
					outFile << atomData[n].atomSum << "\t";
					for (int i = 0; i < AtomData::PARAM_NUM; i++) {
						outFile << atomData[n].param[i] << "\t";
					}
					int moveValue = 0;
					for (int i = 0; i < AtomData::PARAM_NUM; i++) {
						moveValue *= 10;
						moveValue += atomData[n].move[i];
					}
					outFile << moveValue << "\t";
					for (int i = 0; i < AtomData::PARAM_NUM; i++) {
						string connectAtom = "";
						string connectPar = "";
						string connectValue = "";
						for (int j = 0; j < tmpNum[n][i]; j++) {
							connectAtom = connectAtom + "," + to_string(tmpIndex[n][i][j]);
							connectPar = connectPar + "," + to_string(tmpParam[n][i][j]);
							connectValue = connectValue + "," + to_string(tmpRatio[n][i][j]);
						}
						if (connectAtom.length() > 0) {
							outFile << connectAtom.substr(1) << "\t";
						}
						else {
							outFile << "\t";
						}
						if (connectPar.length() > 0) {
							outFile << connectPar.substr(1) << "\t";
						}
						else {
							outFile << "\t";
						}
						if (connectValue.length() > 0) {
							outFile << connectValue.substr(1) << "\t";
						}
						else {
							outFile << "\t";
						}

					}
					for (int i = 0; i < AtomData::PARAM_NUM; i++) {
						outFile << atomData[n].prospectMin[i] << "\t" << atomData[n].prospectMax[i] << "\t";
					}
					outFile << std::endl;

				}
			}
			catch (exception e) {
				outError(e, __LINE__);
			}
			catch (...) {
				runtime_error e("Unknown Error");
				outError(e, __LINE__);
			}
		}
		complex<double> calcF(const int index, const  double h, const double k, const double l, const double q, const complex<double>f) const {
			return atomData[index].calcF(h, k, l, q, f);

		}
		double movePar(Replica* replica, const int atomIndex, const int parIndex, const double dx);
	};
	class ExpData {
	public:
		static const int MAX_DAT = 3000;
		double h[MAX_DAT] = { 0 };
		double k[MAX_DAT] = { 0 };
		double l[MAX_DAT] = { 0 };
		double Iexp[MAX_DAT] = { 0 };
		double Istat[MAX_DAT] = { 0 };
		double q[MAX_DAT] = { 0 };
		complex<double> f[MAX_DAT][FilmData::MAX_ATOM_NUM] = { 0 };
		int dataNum = 0;
		void storeF(FilmData, BulkData, F0Data&);
		ExpData(const string inputFileName) {
			ifstream inputFile(inputFileName);
			try {
				if (!inputFile) {
					throw runtime_error(string("FileNotFoundException(") + inputFileName + string(")"));
				}
				string str;
				while (getline(inputFile, str)) {
					if (dataNum > MAX_DAT) {
						throw runtime_error(string("too many intensity datapoints->" + to_string(dataNum + 1) + string("MAX_DAT_NUM")));
					}
					try {
						vector<string>strList = splitStr(str, '\t');
						if (strList.size() < 4) {
							outLog("Insufficient colum number in the intensity data file. The software ignores the line. (" + str + ")");
							continue;
						}
						strList.resize(6);
						h[dataNum] = stod(strList[0]);
						k[dataNum] = stod(strList[1]);
						l[dataNum] = stod(strList[2]);
						Iexp[dataNum] = stod(strList[3]);
						try {
							Istat[dataNum] = stod(strList[5]);
						}
						catch (invalid_argument& e) {
							outLog(e.what());
							outLog("Because of the error, Istat is set to default value (0)");
							Istat[dataNum] = 0;
						}
					}
					catch (exception& e) {
						outLog(string("Check the line:") + str);
						outError(e, __LINE__);
						throw e;
					}
					catch (...) {
						runtime_error e("Unknown error");
						outLog(string("Check the line:") + str);
						outError(e, __LINE__);
						throw e;
					}
					dataNum++;
				}
			}
			catch (exception& e) {
				outError(e, __LINE__);
				throw e;
			}
			inputFile.close();
		}
	};
	class BulkData {
	public:
		static const int MAX_ATOM_NUM = 30;
		int atomNum = 0;
		AtomData atoms[MAX_ATOM_NUM];
		double a = 0, b = 0, c = 0, alpha = 0, beta = 0, gamma = 0;// angle is written in the unit of radian.  角度はラジアンで統一
		double S11 = 0, S22 = 0, S33 = 0, S12 = 0, S23 = 0, S31 = 0, V = 0;
		/*
		バルクファイルの入力処理
		バルクファイルはタブ区切りで書式は
		１行目:a b c alpha beta gamma
		２行目以降:通し番号 ダミー 原子番号 x y z B occ
		Read out the bulk file.
		1st line: a b c alpha beta gamma
		2nd line and later: index dummy atomic_number x y z B occ
		*/
		BulkData() {};
		BulkData(const string inputBulkFileName) {
			ifstream inputBulkFile(inputBulkFileName);
			try {
				if (!inputBulkFile) {
					throw runtime_error(string("FileNotFoundException(") + inputBulkFileName + string(")"));
				}

				string str;
				getline(inputBulkFile, str);
				vector<string>bulk1 = splitStr(str, '\t');
				if (bulk1.size() < 6) {
					runtime_error e("Bulk file's first line has error");
					throw e;
				}
				a = stod(bulk1[0]);
				b = stod(bulk1[1]);
				c = stod(bulk1[2]);
				alpha = stod(bulk1[3]) * PI / 180;
				beta = stod(bulk1[4]) * PI / 180;
				gamma = stod(bulk1[5]) * PI / 180;
				V = a * b * c * sqrt(1 - pow(cos(alpha), 2) - pow(cos(beta), 2) - pow(cos(gamma), 2) + 2 * cos(alpha) * cos(beta) * cos(gamma));
				S11 = b * b * c * c * pow(sin(alpha), 2);
				S22 = a * a * c * c * pow(sin(beta), 2);
				S33 = b * b * a * a * pow(sin(gamma), 2);
				S12 = a * b * c * c * (cos(alpha) * cos(beta) - cos(gamma));
				S23 = a * a * b * c * (cos(gamma) * cos(beta) - cos(alpha));
				S31 = a * b * b * c * (cos(alpha) * cos(gamma) - cos(beta));

				while (getline(inputBulkFile, str)) {
					try {
						if (atomNum >= MAX_ATOM_NUM) {
							throw runtime_error(string("Too many bulk atoms.->" + to_string(atomNum + 1) + string("MAX_ATOM_NUM")));
						}
						vector<string>strList = splitStr(str, '\t');
						if (strList.size() < 3) {
							outLog("Insufficient column number in the bulk file. The software ignores the line. (" + str + ")");
							continue;
						}
						strList.resize(8);
						atoms[atomNum].atomicNumber = stoi(strList[2]);
						atoms[atomNum].param[0] = 1;
						for (long i = 0; i < 5; i++) {
							atoms[atomNum].param[i + 1] = stod(strList[3 + i]);
						}
						atoms[atomNum].param[6] = 1;
					}
					catch (exception& e) {
						outLog(string("Check the line:") + str);
						outError(e, __LINE__);
						throw e;
					}
					catch (...) {
						runtime_error e("Unknown error");
						outLog(string("Check the line:") + str);
						outError(e, __LINE__);
						throw e;
					}
					atomNum++;
				}
				inputBulkFile.close();
			}
			catch (exception& e) {
				outError(e, __LINE__);
				throw e;
			}
			catch (...) {
				runtime_error e("Unknown error");
				outError(e, __LINE__);
				throw e;

			}

		}
		double calcQ(const double h, const double k, const double l) const {
			try {
				return sqrt((S11 * h * h + S22 * k * k + S33 * l * l + 2 * S12 * h * k + 2 * S23 * k * l + 2 * S31 * h * l) / (V * V)) / 2;
			}
			catch (exception& e) {
				outError(e, __LINE__);
				throw e;
			}
			catch (...) {
				runtime_error e("Unknown error");
				outError(e, __LINE__);
				throw e;

			}
		}


		//ここでのBはBulk本来の温度因子に掛け合わせるfactor
		// The parameter Bfactor here is the multiplier for the bulk atomic displacement parameter.
		complex<double>calcF(const double q, const double h, const double k, const double l, const  F0Data f0Data, const double Bfactor) const {
			try {
				complex<double>F = 0;
				for (int i = 0; i < atomNum; i++) {
					F += f0Data.f(atoms[i].atomicNumber, q) * atoms[i].param[AtomData::toccf] * atoms[i].param[AtomData::foccf]
						* exp(2.0 * PI * I * (h * atoms[i].param[AtomData::x] + k * atoms[i].param[AtomData::y] + l * (atoms[i].param[AtomData::z] - 1.0)))
						* exp(-atoms[i].param[AtomData::B] * Bfactor * q * q) / (1.0 - exp(-2.0 * PI * I * l - MU)) * exp(MU * (atoms[i].param[AtomData::z] - 1.0));
				}

				return F;
			}
			catch (exception& e) {
				outError(e, __LINE__);
				throw e;
			}
			catch (...) {
				runtime_error e("Unknown error");
				outError(e, __LINE__);
				throw e;

			}

		}
	};

	/*
	レプリカモンテカルロ法の一つのレプリカに対応する
	*/
	class Replica {
	public:
		BulkData bulk;
		ExpData* expData;
		F0Data* f0Data;
		FilmData film;
		complex<double>F[ExpData::MAX_DAT] = { 0 };
		complex<double>Ff[FilmData::MAX_ATOM_NUM][ExpData::MAX_DAT] = { 0 };
		double E = 0;

		double data[LOOP_NUM][FilmData::MAX_ATOM_NUM][AtomData::PARAM_NUM] = { 0 };
		int Tindex;
		double revT = 0;// inverse temperature 逆温度
		void outFile(const string outName) const {
			film.outPosFile(outName + "-" + to_string(Tindex) + ".pos");
			outDatFile(outName + "-" + to_string(Tindex) + ".dat");
		}
		void outDatFile(const string outName) const {
			try {
				ofstream outFile(outName);
				for (int n = 0; n < expData->dataNum; n++) {
					outFile << expData->h[n] << "\t" << expData->k[n] << "\t" << expData->l[n] << "\t";
					outFile << expData->Iexp[n] << "\t";
					outFile << film.scale * (F[n].real() * F[n].real() + F[n].imag() * F[n].imag()) << "\t";
					outFile << expData->Istat[n] << std::endl;
				}
			}
			catch (exception& e) {
				outError(e, __LINE__);
				throw e;
			}
			catch (...) {
				runtime_error e("Unknown error");
				outError(e, __LINE__);
				throw e;
			}
		}
		void calcF() {

			try {
				expData->storeF(film, bulk, *f0Data);
				for (int n = 0; n < expData->dataNum; n++) {
					F[n] = bulk.calcF(expData->q[n], expData->h[n], expData->k[n], expData->l[n], *f0Data, film.bulkB);
					for (int i = 0; i < film.atomNum; i++) {
						Ff[i][n] = film.calcF(i, expData->h[n], expData->k[n], expData->l[n], expData->q[n], expData->f[n][i]);
						F[n] += Ff[i][n];
					}
				}
			}
			catch (exception& e) {
				outError(e, __LINE__);
				throw e;
			}
			catch (...) {
				runtime_error e("Unknown error");
				outError(e, __LINE__);
				throw e;
			}

		}
		void renewF(int index) {
			try {
				for (int n = 0; n < expData->dataNum; n++) {
					F[n] -= Ff[index][n];
					Ff[index][n] = film.calcF(index, expData->h[n], expData->k[n], expData->l[n], expData->q[n], expData->f[n][index]);
					F[n] += Ff[index][n];
				}
			}
			catch (exception& e) {
				outError(e, __LINE__);
				throw e;
			}
			catch (...) {
				runtime_error e("Unknown error");
				outError(e, __LINE__);
				throw e;
			}

		}

		void renewF_B(double Bfactor) {
			try {
				for (int n = 0; n < expData->dataNum; n++) {
					F[n] = bulk.calcF(expData->q[n], expData->h[n], expData->k[n], expData->l[n], *f0Data, Bfactor);
					for (int i = 0; i < film.atomNum; i++) {
						F[n] += Ff[i][n];
					}
				}
			}
			catch (exception& e) {
				outError(e, __LINE__);
				throw e;
			}
			catch (...) {
				runtime_error e("Unknown error");
				outError(e, __LINE__);
				throw e;
			}

		}
		double calcE() const {
			try {
				double Esum = 0;
				for (int n = 0; n < expData->dataNum; n++) {
					Esum += pow(expData->Iexp[n] - (film.scale * abs(F[n]) * abs(F[n])), 2)
						/ 2 / (pow(film.noiseScale * (film.scale * abs(F[n]) * abs(F[n])), 2) + pow(film.noiseValue + expData->Istat[n], 2))
						+ 0.5 * log(2 * PI) + 0.5 * log(pow(film.noiseScale * (film.scale * abs(F[n]) * abs(F[n])), 2) + pow(film.noiseValue + expData->Istat[n], 2));
				}
				return Esum / expData->dataNum;
			}
			catch (exception& e) {
				outError(e, __LINE__);
				throw e;
			}
			catch (...) {
				runtime_error e("Unknown error");
				outError(e, __LINE__);
				throw e;
			}
		}
		void MCMC(int loop) {
			if (revT < 0) {
				return;
			}
			E = 0;
			try {
				double T = revT;
				for (int l = 0; l < loop; l++) {
					int fitIndex = 0;
					int atomIndex = 0;
					double Eold = calcE();
					double Enew = 0;
					double Pold = 0;
					double Pnew = 0;
					while (true) {
						atomIndex = rand() % (film.atomNum + 1);
						if (atomIndex == film.atomNum) {
							fitIndex = rand() % 4 + AtomData::PARAM_NUM;
							if (fitIndex == AtomData::PARAM_NUM && film.scaleFit) {
								break;
							}
							if (fitIndex == AtomData::PARAM_NUM + 1 && film.noiseScaleFit) {
								break;
							}
							if (fitIndex == AtomData::PARAM_NUM + 2 && film.noiseValueFit) {
								break;
							}
							if (fitIndex == AtomData::PARAM_NUM + 3 && film.bulkBfit) {
								break;
							}
						}
						else {
							fitIndex = rand() % AtomData::PARAM_NUM;
							if (film.atomData[atomIndex].move[fitIndex]) {
								break;
							}
						}
					}

					if (fitIndex < AtomData::PARAM_NUM) {
						double dr = DELTA * gauss_rand(true);
						if (fitIndex >= 4) {
							dr = dr * 5;
						}
						dr = film.movePar(this, atomIndex, fitIndex, dr);
						Enew = calcE();
						Pnew = 0;
						double accept = exp((-Enew + Eold) * T - Pnew + Pold);
						if (accept < gauss_rand(false)) {
							film.movePar(this, atomIndex, fitIndex, -dr);
						}
					}

					if (fitIndex == AtomData::PARAM_NUM) {
						double oldScale = film.scale;
						film.scale = film.scale * (1 + 0.01 * gauss_rand(true));
						Enew = calcE();
						double accept = exp((-Enew + Eold) * T);
						if (accept < gauss_rand(false)) {
							film.scale = oldScale;
						}
					}
					if (fitIndex == AtomData::PARAM_NUM + 1) {
						double oldNoiseScale = film.noiseScale;
						film.noiseScale = film.noiseScale * (1 + 0.05 * gauss_rand(true));
						Enew = calcE();

						double accept = exp((-Enew + Eold) * T);
						if (accept < gauss_rand(false)) {
							film.noiseScale = oldNoiseScale;
						}
					}
					if (fitIndex == AtomData::PARAM_NUM + 2) {
						double oldNoiseValue = film.noiseValue;
						film.noiseValue = film.noiseValue * (1 + 0.05 * gauss_rand(true));
						Enew = calcE();

						double accept = exp((-Enew + Eold) * T);
						if (accept < gauss_rand(false)) {
							film.noiseValue = oldNoiseValue;
						}
					}
					if (fitIndex == AtomData::PARAM_NUM + 3) {
						double oldB = film.bulkB;
						double dB = gauss_rand(true);
						if (oldB + dB < 0) {

						}
						else {
							renewF_B(oldB + dB);
							Enew = calcE();

							double accept = exp((-Enew + Eold) * T);
							if (accept < gauss_rand(false)) {
								renewF_B(oldB);
							}
							else {
								film.bulkB = oldB + dB;
							}


						}

					}
					for (int n = 0; n < film.atomNum; n++) {
						for (int i = 0; i < AtomData::PARAM_NUM; i++) {
							data[l][n][i] = film.atomData[n].param[i];
						}
					}
					E += Eold;
				}

			}
			catch (exception& e) {
				outError(e, __LINE__);
				throw e;
			}
			catch (...) {
				runtime_error e("Unknown error");
				outError(e, __LINE__);
				throw e;
			}
		}
		static int repIndexFromTindex(Replica rep[], int Tindex) {
			for (int T = 0; T < REPLICA_NUM; T++) {
				if (rep[T].Tindex == Tindex) {
					return T;
				}
			}
			return -1;
		}


	};
	void ExpData::storeF(FilmData film, BulkData bulk, F0Data& f0) {
		try {
			for (int data = 0; data < dataNum; data++) {
				q[data] = bulk.calcQ(h[data], k[data], l[data]);

				for (int atom = 0; atom < film.atomNum; atom++) {

					f[data][atom] = f0.f(film.atomData[atom].atomicNumber, q[data]);
				}
			}
		}
		catch (exception& e) {
			outError(e, __LINE__);
			throw e;
		}
		catch (...) {
			runtime_error e("Unknown error");
			outError(e, __LINE__);
			throw e;
		}
	}
	double FilmData::movePar(Replica* replica, int atomIndex, int parIndex, double dx) {
		try {
			if (atomData[atomIndex].param[parIndex] + dx < atomData[atomIndex].prospectMin[parIndex]) {
				dx = 0;
			}
			if (atomData[atomIndex].param[parIndex] + dx > atomData[atomIndex].prospectMax[parIndex]) {
				dx = 0;
			}

			atomData[atomIndex].param[parIndex] += dx;
			replica->renewF(atomIndex);
			for (int i = 0; i < atomData[atomIndex].connectNum[parIndex]; i++) {
				atomData[atomData[atomIndex].connectAtomIndex[parIndex][i]].param[atomData[atomIndex].connectParIndex[parIndex][i]] += atomData[atomIndex].connectRatio[parIndex][i] * dx;
				replica->renewF(atomData[atomIndex].connectAtomIndex[parIndex][i]);
			}
			return dx;
		}
		catch (exception& e) {
			outError(e, __LINE__);
			throw e;
		}
		catch (...) {
			runtime_error e("Unknown error");
			outError(e, __LINE__);
			throw e;
		}


	}

	/*
	最初に呼ばれるmain関数
	それぞれの引数の数に応じたmain関数を呼び出す。
	*/
	int main(int argc, char* argv[]) {
		try {
			logFile = ofstream("log.txt", std::ios::app);
			std::cout << "CTR-structure, Tohoku Univ." << std::endl;

			if (argc < 6) {
				std::cout << "Insufficient arguments" << std::endl;
				std::cout << "Sample arguments : CTR-structure-3d.exe inputBulkFileName inputFilmFileName inputExpFileName outPutFileName iter maxBeta minBeta";
				throw runtime_error("Insufficient arguments");
			}
			string inputBulkFileName = string(argv[1]);
			string inputFilmFileName = string(argv[2]);
			string inputExpFileName = string(argv[3]);
			string outPutFileName = string(argv[4]);
			int iter = atoi(argv[5]);
			double minBeta = 0.5;
			double maxBeta = 500;
			if (argc > 7) {
				minBeta = atof(argv[7]);
				maxBeta = atof(argv[6]);
			}
			else {
				outLog("No beta input. Software use maxBeta=500 minBeta=0.5");
			}



			outLog("");
			outLog("<start CTR-structure>");
			outLog("start time\t" + getTime());
			outLog("output file(folder) name\t" + outPutFileName);
			outLog("bulk's filename\t" + inputBulkFileName);
			outLog("film filename\t" + inputFilmFileName);
			outLog("intensity filename\t" + inputExpFileName);
			outLog("DELTA\t" + to_string(DELTA));
			outLog("LOOP_NUM\t" + to_string(LOOP_NUM));
			outLog("MU\t" + to_string(MU));


			F0Data f0Data("f0data.dat");
			f0Data.anomalousLoad("anomalous.dat");

			BulkData bulkData(inputBulkFileName);
			FilmData* filmData;
			filmData = new FilmData(inputFilmFileName);
			ExpData* expData;
			expData = new ExpData(inputExpFileName);


			if (iter == 0) {
				Replica* replica;
				replica = new Replica;
				replica->film = *filmData;
				replica->bulk = bulkData;
				replica->expData = expData;
				replica->f0Data = &f0Data;
				replica->revT = 1;
				replica->Tindex = 0;
				replica->calcF();
				replica->film.outPosFile(outPutFileName + ".pos");
				replica->outDatFile(outPutFileName + ".dat");
				return 0;

			}

			if (_mkdir(outPutFileName.c_str()) != 0) {
				runtime_error e("output folder already exists.");
				throw e;
			}

			if (_chdir(outPutFileName.c_str()) != 0) {
				runtime_error e("can't open output folder");
				throw e;
			}
			Replica* replica;
			replica = new Replica[REPLICA_NUM];
			double reverseTlist[REPLICA_NUM];

			ofstream ElogFile(outPutFileName + ".log");
			ofstream swapFile(outPutFileName + ".swap");

			for (int T = 0; T < REPLICA_NUM; T++) {
				replica[T].film = *filmData;
				replica[T].bulk = bulkData;
				replica[T].expData = expData;
				replica[T].f0Data = &f0Data;
				reverseTlist[T] = std::exp(T / (double)(REPLICA_NUM - 1) * (std::log(minBeta / maxBeta)) + std::log(maxBeta));
				replica[T].revT = reverseTlist[T];
				replica[T].Tindex = T;
				replica[T].calcF();
			}
			replica[0].film.outPosFile("startFile.pos");
			replica[0].outDatFile("startFile.dat");

			ofstream outPar[REPLICA_NUM];
			for (int T = 0; T < REPLICA_NUM; T += 8) {
				outPar[T].open(outPutFileName + "-" + to_string(T) + ".par");
			}
			outLog("success inializing and loading files.");
			outLog("start loop calculation");
			for (int i = 0; i < iter / LOOP_NUM; i++) {
#pragma omp parallel 
				{
#pragma omp for
					for (int T = 0; T < REPLICA_NUM; T++) {
						replica[T].MCMC(LOOP_NUM);
					}
				}
				try {
					// Exchange replicas レプリカ交換
					for (int T = REPLICA_NUM - 1; T > 0; T = T - 2) {
						int R1 = Replica::repIndexFromTindex(replica, T);
						int R2 = Replica::repIndexFromTindex(replica, T - 1);
						if (replica[R1].revT < 0) {
							continue;
						}
						double exchangeP = std::exp((replica[R1].revT - replica[R2].revT) * (replica[R1].calcE() - replica[R2].calcE()));
						if (exchangeP > gauss_rand(false)) {
							replica[R1].Tindex = T - 1;
							replica[R2].Tindex = T;
							double tmpRevT = replica[R1].revT;
							replica[R1].revT = replica[R2].revT;
							replica[R2].revT = tmpRevT;

						}
					}
					for (int T = REPLICA_NUM - 2; T > 0; T = T - 2) {
						int R1 = Replica::repIndexFromTindex(replica, T);
						int R2 = Replica::repIndexFromTindex(replica, T - 1);
						if (replica[R1].revT < 0) {
							continue;
						}
						double exchangeP = std::exp((replica[R1].revT - replica[R2].revT) * (replica[R1].calcE() - replica[R2].calcE()));
						if (exchangeP > gauss_rand(false)) {
							replica[R1].Tindex = T - 1;
							replica[R2].Tindex = T;
							double tmpRevT = replica[R1].revT;
							replica[R1].revT = replica[R2].revT;
							replica[R2].revT = tmpRevT;
						}
					}
				}
				catch (exception& e) {
					outError(e, __LINE__);
					throw e;
				}
				catch (...) {
					runtime_error e("Unknown error");
					outError(e, __LINE__);
					throw e;
				}
				try {
					for (int T = 0; T < REPLICA_NUM; T += 8) {
						int replicaIndex = Replica::repIndexFromTindex(replica, T);
						for (int l = 0; l < LOOP_NUM; l += RECORD_NUM) {
							for (int n = 0; n < filmData->atomNum; n++) {
								for (int k = 0; k < AtomData::PARAM_NUM; k++) {
									if (filmData->atomData[n].move[k]) {
										outPar[T] << replica[replicaIndex].data[l][n][k] << "\t";
									}
								}
							}
							outPar[T] << std::endl;
						}
					}
					for (int T = 0; T < REPLICA_NUM; T++) {
						ElogFile << replica[Replica::repIndexFromTindex(replica, T)].E / LOOP_NUM << "\t";
						swapFile << replica[T].Tindex << "\t";
					}
					swapFile << std::endl;

					ElogFile << std::endl;
				}
				catch (exception& e) {
					outError(e, __LINE__);
					throw e;
				}
				catch (...) {
					runtime_error e("Unknown error");
					outError(e, __LINE__);
					throw e;
				}


			}


			for (int T = 0; T < REPLICA_NUM; T++) {
				replica[T].outFile(outPutFileName);
			}
			outLog("finished time\t" + getTime());
			outLog("<success calculation>");
		}
		catch (exception& e) {
			outError(e, __LINE__);
			return -1;
		}
		catch (...) {
			runtime_error e("Unknown error");
			outError(e, __LINE__);
			return -1;
		}
		return 0;

	}

}
//デフォルト実行関数
int main(int argc, char* argv[]) {
	CTRstructure::main(argc, argv);
}