#include "BspMap.h"
#include "Utils.h"
#include <assert.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace KernelBridgeNS;
using namespace SweepNS;
using namespace std;

template <class PTypeD, class PTypeR>
BspMap<PTypeD, PTypeR>::BspMap(const string & fileName)
{
	ifstream
		fileStream(fileName);
	string fileBuffer, stringToken;
	SInt
		bufferInd = 0;

	if (!fileStream.is_open()) {
		printf("Unable to open file: %s\n", fileName.c_str());
		return;
	}

/*
	while ((token = getToken(fileStream, fileBuffer, bufferInd, stringToken)) != TOKEN_EOF && token != TOKEN_OPEN_PAREN);

	if (getToken(fileStream, fileBuffer, bufferInd, stringToken) != IP_TOKEN_OBJECT) {
		printf("OBJECT key words expected\n");
		return;
	}

	if (getToken(fileStream, fileBuffer, bufferInd, stringToken) != TOKEN_OTHER) {
		printf("Name expected\n");
		return;
	}

	if (getToken(fileStream, fileBuffer, bufferInd, stringToken) != TOKEN_OPEN_PAREN) {
		printf("\"[\" expected\n");
		return;
	}

	if (getToken(fileStream, fileBuffer, bufferInd, stringToken) != TOKEN_MULTIVAR ||
		getToken(fileStream, fileBuffer, bufferInd, stringToken) != TOKEN_BSPLINE) {
		printf("MVAR BSPLINE key words expected\n");
		return;
	}

	if ((token = getToken(fileStream, fileBuffer, bufferInd, stringToken)) != TOKEN_OTHER ||
		sscanf_s(stringToken.c_str(), "%d", &_domainDim) != 1) {
		printf("BSPLINE's dimension expected\n");
		return;
	}

	_lengths = vector<SInt>(_domainDim);
	for (SInt i = 0; i < _domainDim; i++) {
		if ((token = getToken(fileStream, fileBuffer, bufferInd, stringToken)) != TOKEN_OTHER ||
			sscanf_s(stringToken.c_str(), "%d", &_lengths[i]) != 1) {
			printf("BSPLINE's lengths of mesh expected");
			return;
		}
	}

	_orders = vector<SInt>(_domainDim);
	for (SInt i = 0; i < _domainDim; i++) {
		if ((token = getToken(fileStream, fileBuffer, bufferInd, stringToken)) != TOKEN_OTHER ||
			sscanf_s(stringToken.c_str(), "%d", &_orders[i]) != 1) {
			printf("BSPLINE's orders expected\n");
			return;
		}
	}

	_subSpaces = vector<SInt>(_domainDim + 1);
	_subSpaces[0] = 1;
	for (SInt i = 0; i < _domainDim; i++) {
		assert(_orders[i] <= _lengths[i]);
		_subSpaces[i + 1] = _subSpaces[i] * _lengths[i];
	}

	token = getToken(fileStream, fileBuffer, bufferInd, stringToken);
	if (token < TOKEN_E1 || token > TOKEN_E9 ||
		stringToken.length() != 2 ||
		stringToken[0] != 'E' ||
		!isdigit(stringToken[1]) ||
		atoi(&stringToken[1]) > 9) {
		printf("BSPLINE Point type expected\n");
		return;
	}

	_rangeDim = atoi(&stringToken[1]);

	assert(_rangeDim == PTypeR::dim());

	_knotVectors = vector<vector<SReal>>(_domainDim);
	for (SInt k = 0; k < _domainDim; k++) {
		SInt
			len = _orders[k] + _lengths[k];

		_knotVectors[k] = vector<SReal>(len);
		if ((token = getToken(fileStream, fileBuffer, bufferInd, stringToken)) != TOKEN_OPEN_PAREN) {
			printf("\"[\" expected\n");
			return;
		}
		if ((token = getToken(fileStream, fileBuffer, bufferInd, stringToken)) != TOKEN_KV) {
			printf("KV/KVP expected\n");
			return;
		}

		for (SInt i = 0; i < len; i++) {
			if ((token = getToken(fileStream, fileBuffer, bufferInd, stringToken)) != TOKEN_OTHER ||
				!getRealNumber(stringToken, _knotVectors[k][i])) {
				printf("Numeric data expected\n");
				return;
			}
		}

		if ((token = getToken(fileStream, fileBuffer, bufferInd, stringToken)) != TOKEN_CLOSE_PAREN) {
			printf("\"]\" expected\n");
			return;
		}
	}

	_points = vector<PTypeR>(_subSpaces[_domainDim]);
	for (SInt i = 0; i < _subSpaces[_domainDim]; i++) {
		if ((token = getToken(fileStream, fileBuffer, bufferInd, stringToken)) != TOKEN_OPEN_PAREN) {
			printf("\"[\" expected\n");
			return;
		}
		_points[i] = PTypeR();
		for (SInt j = 0; j < _rangeDim; j++) {
			double R;

			if ((token = getToken(fileStream, fileBuffer, bufferInd, stringToken)) != TOKEN_OTHER ||
				!getRealNumber(stringToken, R)) {
				printf("Numeric data expected\n");
				return;
			}
			_points[i].setCoord(j, R);
		}
		if ((token = getToken(fileStream, fileBuffer, bufferInd, stringToken)) != TOKEN_CLOSE_PAREN) {
			printf("\"]\" expected\n");
			return;
		}
	}

	if ((token = getToken(fileStream, fileBuffer, bufferInd, stringToken)) != TOKEN_CLOSE_PAREN) {
		printf("\"]\" expected\n");
		return;
	}
*/
		
	readFromStream(fileStream);

	fileStream.close();
}

template <class PTypeD, class PTypeR>
BspMap<PTypeD, PTypeR>::BspMap(ifstream & inFS)
{
	assert(inFS.is_open());		
	readFromStream(inFS);
}


template <class PTypeD, class PTypeR>
void BspMap<PTypeD, PTypeR>::readFromStream(ifstream & inFS)
{
	string fileBuffer, stringToken;
	SInt
		bufferInd = 0;
	ITDTokenType token;

	assert(inFS.is_open());

	while ((token = getToken(inFS, fileBuffer, bufferInd, stringToken)) != TOKEN_EOF && token != TOKEN_OPEN_PAREN);

	if (getToken(inFS, fileBuffer, bufferInd, stringToken) != IP_TOKEN_OBJECT) {
		printf("OBJECT key words expected\n");
		return;
	}

	if (getToken(inFS, fileBuffer, bufferInd, stringToken) != TOKEN_OTHER) {
		printf("Name expected\n");
		return;
	}

	if (getToken(inFS, fileBuffer, bufferInd, stringToken) != TOKEN_OPEN_PAREN) {
		printf("\"[\" expected\n");
		return;
	}

	if (getToken(inFS, fileBuffer, bufferInd, stringToken) != TOKEN_MULTIVAR ||
		getToken(inFS, fileBuffer, bufferInd, stringToken) != TOKEN_BSPLINE) {
		printf("MVAR BSPLINE key words expected\n");
		return;
	}

	if ((token = getToken(inFS, fileBuffer, bufferInd, stringToken)) != TOKEN_OTHER ||
		sscanf_s(stringToken.c_str(), "%d", &_domainDim) != 1) {
		printf("BSPLINE's dimension expected\n");
		return;
	}

	_lengths = vector<SInt>(_domainDim);
	for (SInt i = 0; i < _domainDim; i++) {
		if ((token = getToken(inFS, fileBuffer, bufferInd, stringToken)) != TOKEN_OTHER ||
			sscanf_s(stringToken.c_str(), "%d", &_lengths[i]) != 1) {
			printf("BSPLINE's lengths of mesh expected");
			return;
		}
	}

	_orders = vector<SInt>(_domainDim);
	for (SInt i = 0; i < _domainDim; i++) {
		if ((token = getToken(inFS, fileBuffer, bufferInd, stringToken)) != TOKEN_OTHER ||
			sscanf_s(stringToken.c_str(), "%d", &_orders[i]) != 1) {
			printf("BSPLINE's orders expected\n");
			return;
		}
	}

	_subSpaces = vector<SInt>(_domainDim + 1);
	_subSpaces[0] = 1;
	for (SInt i = 0; i < _domainDim; i++) {
		assert(_orders[i] <= _lengths[i]);
		_subSpaces[i + 1] = _subSpaces[i] * _lengths[i];
	}

	token = getToken(inFS, fileBuffer, bufferInd, stringToken);
	if (token < TOKEN_E1 || token > TOKEN_E9 ||
		stringToken.length() != 2 ||
		stringToken[0] != 'E' ||
		!isdigit(stringToken[1]) ||
		atoi(&stringToken[1]) > 9) {
		printf("BSPLINE Point type expected\n");
		return;
	}

	_rangeDim = atoi(&stringToken[1]);

	assert(_rangeDim == PTypeR::dim());

	_knotVectors = vector<vector<SReal>>(_domainDim);
	for (SInt k = 0; k < _domainDim; k++) {
		SInt
			len = _orders[k] + _lengths[k];

		_knotVectors[k] = vector<SReal>(len);
		if ((token = getToken(inFS, fileBuffer, bufferInd, stringToken)) != TOKEN_OPEN_PAREN) {
			printf("\"[\" expected\n");
			return;
		}
		if ((token = getToken(inFS, fileBuffer, bufferInd, stringToken)) != TOKEN_KV) {
			printf("KV/KVP expected\n");
			return;
		}

		for (SInt i = 0; i < len; i++) {
			if ((token = getToken(inFS, fileBuffer, bufferInd, stringToken)) != TOKEN_OTHER ||
				!getRealNumber(stringToken, _knotVectors[k][i])) {
				printf("Numeric data expected\n");
				return;
			}
		}

		if ((token = getToken(inFS, fileBuffer, bufferInd, stringToken)) != TOKEN_CLOSE_PAREN) {
			printf("\"]\" expected\n");
			return;
		}
	}

	_points = vector<PTypeR>(_subSpaces[_domainDim]);
	for (SInt i = 0; i < _subSpaces[_domainDim]; i++) {
		if ((token = getToken(inFS, fileBuffer, bufferInd, stringToken)) != TOKEN_OPEN_PAREN) {
			printf("\"[\" expected\n");
			return;
		}
		_points[i] = PTypeR();
		for (SInt j = 0; j < _rangeDim; j++) {
			double R;

			if ((token = getToken(inFS, fileBuffer, bufferInd, stringToken)) != TOKEN_OTHER ||
				!getRealNumber(stringToken, R)) {
				printf("Numeric data expected\n");
				return;
			}
			_points[i].setCoord(j, R);
		}
		if ((token = getToken(inFS, fileBuffer, bufferInd, stringToken)) != TOKEN_CLOSE_PAREN) {
			printf("\"]\" expected\n");
			return;
		}
	}

	if ((token = getToken(inFS, fileBuffer, bufferInd, stringToken)) != TOKEN_CLOSE_PAREN) {
		printf("\"]\" expected\n");
		return;
	}

	if ((token = getToken(inFS, fileBuffer, bufferInd, stringToken)) != TOKEN_CLOSE_PAREN) {
		printf("\"]\" expected\n");
		return;
	}
}


template <class PTypeD, class PTypeR>
int BspMap<PTypeD, PTypeR>::getRealNumber(const string & strNum,
	double & realNum)
{
	unsigned int i;

	/* For those that use ',' instead of '.' as decimal point... */
	if (strchr(strNum.c_str(), ',') != NULL && strchr(strNum.c_str(), '.') == NULL) {
		string
			strNum2 = strNum;

		for (i = 0; i < strNum2.length(); i++) {
			if (strNum2[i] == ',')
				strNum2[i] = '.';
		}
		i = sscanf_s(strNum2.c_str(), "%lf", &realNum);
	}
	else
		i = sscanf_s(strNum.c_str(), "%lf", &realNum);

	return i == 1;
}


template <class PTypeD, class PTypeR>
SBool BspMap<PTypeD, PTypeR>::getStringToken(ifstream & fileStream,
	string & fileBuffer,
	SInt & bufferInd,
	string & stringToken,
	SBool & quoted)
{
	char c;
	SInt
		bufferLen = (SInt)fileBuffer.length();

	quoted = SFalse;
	stringToken.clear();

	if (bufferInd >= bufferLen) {
		if (!std::getline(fileStream, fileBuffer))
			return SFalse;

		bufferLen = (SInt)fileBuffer.length();
		bufferInd = 0;
	}

	/* skip white spaces: */
	while ((bufferInd < bufferLen) &&
		(((c = fileBuffer[bufferInd++]) == ' ') ||
		(c == '\t') ||
			(c == '\n') ||
			(c == '\r') ||
			(c == '#')) &&
			(c != (char)EOF)) {
		if (c == '#') {					    /* Skip comment. */
			while (c != '\r' && c != '\n') {
				c = fileBuffer[bufferInd++];
			}
		}
		//if (c == '\n')
		  //  _IPStream[Handler].LineNum++;		 /* Count the lines. */
	}

	if (c == '[') 		      /* Its a token by	itself so return it. */
		stringToken += c;	      /* Copy the token	into string. */
	else {
		if (c != (char)EOF) {
			if (c == '"') {
				quoted = STrue;
				while (bufferInd < bufferLen && ((c = fileBuffer[bufferInd++]) != '"')) {
					stringToken += c;      /* Copy the quoted string. */
					if (c == '\\') {
						/* Next character is quoted - copy verbatim. */
						stringToken[stringToken.length() - 1] = c = fileBuffer[bufferInd++];
					}
				}
			}
			else {
				do {
					stringToken += c;  /* Copy the token into string. */
				} while ((bufferInd < bufferLen) &&
					((c = fileBuffer[bufferInd++]) != ' ') &&
					(c != '\t') &&
					(c != '\n') &&
					(c != '\r') &&
					(c != (char)EOF));
			}
			//if (!InputEOF(Handler) && c == '\n')
			//IPInputUnGetC(Handler, c);            /* Save for next time. */
		}
	}

	/* The following handles the spacial case were we have XXXX] - we must   */
	/* split it	into two token XXXX and	], _IPUnGetToken(']') & return XXXX: */
	SInt
		len = (SInt)stringToken.length();

	if (!quoted &&
		(len == 0 || (stringToken[--len] == ']' && len > 0))) {
		/* Return CloseParan */
		bufferInd--;	 /* Save next token. */
		stringToken = stringToken.substr(0, len);			/* Set end of string on	"]". */
	}

	return STrue;
}


template <class PTypeD, class PTypeR>
ITDTokenType BspMap<PTypeD, PTypeR>::getToken(ifstream & fileStream,
	string & fileBuffer,
	SInt & bufferInd,
	string & stringToken)
{
	static const vector<ITDTokenType> intTokens = {
	TOKEN_OPEN_PAREN,
	TOKEN_CLOSE_PAREN,

	TOKEN_E1,
	TOKEN_E2,
	TOKEN_E3,
	TOKEN_E4,
	TOKEN_E5,
	TOKEN_E6,
	TOKEN_E7,
	TOKEN_E8,
	TOKEN_E9,

	IP_TOKEN_OBJECT,
	TOKEN_BEZIER,
	TOKEN_BSPLINE,
	TOKEN_PTYPE,
	TOKEN_NUM_PTS,
	TOKEN_ORDER,
	TOKEN_KV,
	TOKEN_MULTIVAR,
	TOKEN_NONE,
	};
	static const vector<string> strTokens = {
	"[",
	"]",
	"E1",
	"E2",
	"E3",
	"E4",
	"E5",
	"E6",
	"E7",
	"E8",
	"E9",

	"OBJECT",
	"BEZIER",
	"BSPLINE",
	"PTYPE",
	"NUMPTS",
	"ORDER",
	"KV",
	"MULTIVAR",
	""
	};
	SInt i;
	SBool quoted;

	if (!getStringToken(fileStream, fileBuffer, bufferInd, stringToken, quoted))
		return TOKEN_EOF;

	if (quoted)
		return TOKEN_QUOTED;

	for (i = 0; i < (SInt)strTokens.size(); i++) {
		if (_stricmp(stringToken.c_str(), strTokens[i].c_str()) == 0)
			return intTokens[i];
	}

	return TOKEN_OTHER;			  /* Must be number or name. */
}


template <class PTypeD, class PTypeR>
void BspMap<PTypeD, PTypeR>::writeToFile(const string & FileName, SBool append)	const
{
	
	SInt 
		flag = append ? (ios::out | ios::app) : (ios::out);
	ofstream outFS(FileName, flag);

	if (!outFS.is_open()) {
		printf("Cannot open %s for writing.\n", FileName.c_str());
		return;
	}
		
	writeToStream(outFS);
	outFS.close();
}

template <class PTypeD, class PTypeR>
void BspMap<PTypeD, PTypeR>::writeToStream(ofstream & outFS)	const
{
	SInt i, j, len;
	char buffer[512];

	sprintf_s(buffer, 512, "[OBJECT %d_%d\n", _id1, _id2);
	outFS << buffer;

	sprintf_s(buffer, 512, "    [MULTIVAR BSPLINE %d  ", _domainDim);
	outFS << buffer;

	for (i = 0; i < _domainDim; i++) {
		sprintf_s(buffer, 512, " %d", _lengths[i]);
		outFS << buffer;
	}
	outFS << "  ";

	for (i = 0; i < _domainDim; i++) {
		sprintf_s(buffer, 512, " %d", _orders[i]);
		outFS << buffer;
	}

	sprintf_s(buffer, 512, " E%d\n", _rangeDim);
	outFS << buffer;

	for (i = 0; i < _domainDim; i++) {
		len = _lengths[i] + _orders[i];
		outFS << "        [KV";
		for (j = 0; j < len; j++) {
			if (j && j % 5 == 0)
				outFS << "\n";

			sprintf_s(buffer, 512, " %-16.14lg", _knotVectors[i][j]);
			string tmpStr(buffer);
			tmpStr.erase(tmpStr.find_last_not_of(" ") + 1);
			outFS << tmpStr;
			//outFS << buffer;
		}
		outFS << "]\n";
	}

	for (i = 0; i < _subSpaces[_domainDim]; i++) {
		outFS << "        [";
		for (j = 0; j < _rangeDim; j++) {
			sprintf_s(buffer, 512, " %-16.14lg", _points[i].coord(j));			
			outFS << buffer;
		}
		outFS << "]\n";
	}

	outFS << "    ]\n";
	outFS << "]\n";
}

template class BspMap<Pnt1D, Pnt1D>;
template class BspMap<Pnt1D, Pnt2D>;
template class BspMap<Pnt1D, Pnt3D>;
template class BspMap<Pnt1D, Pnt4D>;
template class BspMap<Pnt1D, Pnt5D>;
template class BspMap<Pnt1D, Pnt6D>;

template class BspMap<Pnt2D, Pnt1D>;
template class BspMap<Pnt2D, Pnt2D>;
template class BspMap<Pnt2D, Pnt3D>;
template class BspMap<Pnt2D, Pnt4D>;
template class BspMap<Pnt2D, Pnt5D>;
template class BspMap<Pnt2D, Pnt6D>;

template class BspMap<Pnt3D, Pnt1D>;
template class BspMap<Pnt3D, Pnt2D>;
template class BspMap<Pnt3D, Pnt3D>;
template class BspMap<Pnt3D, Pnt4D>;
template class BspMap<Pnt3D, Pnt5D>;
template class BspMap<Pnt3D, Pnt6D>;

template class BspMap<Pnt4D, Pnt1D>;
template class BspMap<Pnt4D, Pnt2D>;
template class BspMap<Pnt4D, Pnt3D>;
template class BspMap<Pnt4D, Pnt4D>;
template class BspMap<Pnt4D, Pnt5D>;
template class BspMap<Pnt4D, Pnt6D>;

template class BspMap<Pnt5D, Pnt1D>;
template class BspMap<Pnt5D, Pnt2D>;
template class BspMap<Pnt5D, Pnt3D>;
template class BspMap<Pnt5D, Pnt4D>;
template class BspMap<Pnt5D, Pnt5D>;
template class BspMap<Pnt5D, Pnt6D>;

template class BspMap<Pnt6D, Pnt1D>;
template class BspMap<Pnt6D, Pnt2D>;
template class BspMap<Pnt6D, Pnt3D>;
template class BspMap<Pnt6D, Pnt4D>;
template class BspMap<Pnt6D, Pnt5D>;
template class BspMap<Pnt6D, Pnt6D>;