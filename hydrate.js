// Version: 0.1
// Author: KMac
// Date: Sept 27, 2021
// Works with https://ipfs.io/ipfs/QmTn16U9YtbMeGkYWtRFZ5XrNA7tiKQYyRSgX4Es3TWouB test data
// Use by calling hydrateLibrary(,,) with the SIPMath Library JSON, Token Name, and number of trials eg 1000

import { jStat } from 'jstat';

export default function hydrateLibrary(LibraryIn, tokenIn, numberOfTrialsRequested) {
  var output = simulateSIP(LibraryIn, tokenIn); // see siplibs/test.json for token names
  let merged = flatten(output);
  return merged.slice(0, numberOfTrialsRequested);
}

const maxTrials = 1000;
function simulateSIP(selfIn, sip) {
  /**
   * Expects library as input and the name of A sip
   * TODO: Add an all option for doing all sips
   * example keyword argument for input: simulateSIP("Variable1", HDR2= {"seed3":0, "seed4":1})
   */

  let randomarray = [];
  let returnValue = [];
  let sipIndex = selfIn.sips.findIndex(item => item.name === sip);

  let aCoeffs = selfIn.sips[sipIndex].arguments.aCoefficients;
  let lowerBound = "";
  let upperBound = "";
  let a = "";
  let functionName = selfIn.sips[sipIndex].function;

  // This need to be rethought out as we are pulling a lot of numbers we don't need but for now :shrugs:
  if (selfIn.sips[sipIndex].ref.source === "copula") {
    randomarray = generateCopula(selfIn, selfIn.sips[sipIndex].ref.copulaLayer); // c1 or c2 etc
  } else if (selfIn.sips[sipIndex].ref.source === "rng") {
    randomarray = prepGenerateRandom(selfIn.sips[sipIndex].ref.name, selfIn.U01);
  }

  try {
    lowerBound = selfIn.sips[sipIndex]["arguments"]["lowerBound"];
  } catch (e) {
    //km console.log("Nothing to see here. Just no lowerBound");
  }

  try {
    upperBound = selfIn.sips[sipIndex]["arguments"]["upperBound"];
  } catch (e) {
    //km console.log("Nothing to see here. Just no upperBound");
  }

  if (functionName === "Metalog_1_0") {
    // TODO: change for loop into a vectorized numpy like function.
    if (randomarray.constructor === Array) {
      for (var i = 0; i < randomarray.length; i++) {
        let ml = metalog(randomarray[i], aCoeffs, lowerBound, upperBound);
        returnValue.push(ml);
      }
    }

    for (let index = 0; index < randomarray.length; ++index) {
      let ml = metalog(randomarray[index], aCoeffs, lowerBound, upperBound);
      returnValue.push(ml);
    }
  }
  return returnValue;
}

function metalog(y, a, bl = "", bu = "") {
  function convert_to_float(a) {
    var floatValue = parseFloat(a);
    return floatValue;
  }

  let vector = [];
  let np_a = a;
  for (let index = 1; index < a.length + 1; ++index) {
    //cant start with 0 coeeffs for each aCoeff
    vector.push(basis(y, index));
  }

  let wrappedVector = [vector];
  let wrappedNp_x = [np_a];

  let wrappedNp_a = wrappedNp_x[0].map(e => [e])
  let mky = multiply(wrappedVector, wrappedNp_a);

  // Unbounded
  if (typeof bl == String && typeof bu == String) {
    return mky;
  }
  if ((typeof bl === "string" || bl instanceof String) && (typeof bu === "string" || bu instanceof String)) {
    return mky;
  }
  // Bounded lower
  else if (typeof bl !== "string" && typeof bu == "string") {
    convert_to_float(bl);
    return bl + Math.exp(mky);
  }
  // Bounded upper
  else if (typeof bl == "string" && typeof bu != "string") {
    convert_to_float(bu);
    return bu - Math.exp(-mky);
  }
  // Bounded
  else if (typeof bl != "string" && typeof bu != "string") {
    return bl + (bu * Math.exp(mky)) / (1 + Math.exp(mky));
  }
}

function generateCopula(selfy, copulaCount) {
  let ret = [];
  let whichCorrelationMatrix = [];

  selfy.U01.copula.forEach((copula) => {
    if (copula.function === "GaussianCopula") {
      // now get the cholesky factors and the global variable
      let specifyCorrelationMatrix = copula.arguments.correlationMatrix.value;
      let copulaArgs = copula.arguments.rng;
      let randomMatrixOfHDRs = [];
      for (let i = 0; i < copulaArgs.length; i++) {
        let val = prepGenerateRandom(copulaArgs[i], selfy.U01); // from U01/RNG
        /* TODO update HDRv2 using { "counter": "PM_Index","entity": 1,"varId": 6187319,"seed3": 0,"seed4": 0} */
        randomMatrixOfHDRs.push(val);
      }

      selfy.globalVariables.forEach((item, index) => {
        if (item["name"] == specifyCorrelationMatrix) {
          whichCorrelationMatrix = index;
        } else {
          let index = -1;
        }
      });

      let thisCorrelationMatrix = selfy.globalVariables[whichCorrelationMatrix].value.matrix;
      let correlationMatrix = convertMx(thisCorrelationMatrix);

      // Find the Cholesky Factors
      let cho = jStat(jStat.cholesky(correlationMatrix));
      //Apply the Cholesky Factors to the randoms
      let col = copula.copulaLayer.findIndex((item) => item === copulaCount);
      let choSubSample = cho[col].slice(0, col + 1);
      let runiRow = [];
      let corrSamples = [];
      let normSinv = [];

      for (let i = 0; i < randomMatrixOfHDRs[0].length; i++) {
        let randomMatrixHRDsSample = [];
        for (let index = 0; index < col + 1; ++index) {
          //each variable upto pos col
          // get first x cols in randuniframe
          randomMatrixHRDsSample[index] = randomMatrixOfHDRs[index];
          runiRow[i] = randomMatrixHRDsSample.map(function (x) {
            return x[i];
          });
        }

        normSinv = runiRow[i].map((sin) => jStat.normal.inv(sin, 0, 1));
        corrSamples[i] = jStat.dot(normSinv, choSubSample); // TODO: Replace jstat 
        corrSamples[i] = jStat.normal.cdf(corrSamples[i], 0, 1); // TODO: Replace jstat 
      }

      ret = corrSamples;
    } else {
      console.log("TypeError The function type for this copula is unsupported.");
    }
  });
  return ret;
}

function prepGenerateRandom(args, selfIn) {
  // from U01/RNG
  let rngArgs = selfIn.rng.findIndex((x) => x.name === args);
  var samples = [];
  const seedPerDist = selfIn.rng[rngArgs].arguments.varId;
  for (var distTrials = 0; distTrials < maxTrials; distTrials++) {
    samples[distTrials] = HDRando(seedPerDist, distTrials);
  }
  return samples;
}

/*
 * hubbardresearch.com for more info. This is a function that generates the random numbers with seeds.
 * TODO update this to use all the seeds from the U01/RNG ie use HRDv2. Move into own package?
 */
function HDRando(seed, PM_Index) {
  const largePrime = 2147483647;
  const million = 1000000;
  const tenMillion = 10000000;

  function mod(n, m) {
    var remain = n % m;
    return Math.floor(remain >= 0 ? remain : remain + m);
  }

  let randi = (
    mod(
      (mod((seed + million) ^ (2 + (seed + million) * (PM_Index + tenMillion)), 99999989) + 1000007) * 
        (mod(
          (PM_Index + tenMillion) ^
            (2 + (PM_Index + tenMillion) *
                mod((seed + million) ^ (2 + (seed + million) * (PM_Index + tenMillion)), 99999989)),
          99999989
        ) + 1000013),
      largePrime) 
      + 0.5) / largePrime;
  return randi;
}

// HELPER FUNCTIONS TODO: Remove need for jstat
function convertMx(correlationMatrix) {
  let variables = [];

  //gotta figure out all of the variables in the matrix
  correlationMatrix.forEach(sipVar => {
    if (variables.includes(sipVar["row"])) {
    } else {
      variables.push(sipVar["row"]);
    }
  });

  let variableCount = variables.length;
  let returnArray = Array(variableCount)
    .fill()
    .map(() => Array(variableCount).fill(0));

  correlationMatrix.forEach(items => {
    let i = items.row;
    let j = items.col;
    let value = items.value;

    i = variables.indexOf(items["row"]);
    j = variables.indexOf(items["col"]);
    returnArray[i][j] = value;
    returnArray[j][i] = value;
  });

  return returnArray;
}

function multiply(a, b) {
  var aNumRows = a.length,
    aNumCols = a[0].length || 0, // if a is a vector
    bNumRows = b.length,
    bNumCols = b[0].length || 0,
    m = new Array(aNumRows); // initialize array of rows

  for (var r = 0; r < aNumRows; ++r) {
    m[r] = new Array(bNumCols); // initialize the current row
    for (var c = 0; c < bNumCols; ++c) {
      m[r][c] = 0; // initialize the current cell
      for (var i = 0; i < aNumCols; ++i) {
        m[r][c] += a[r][i] * b[i][c];
      }
    }
  }

  return m;
}

function basis(y, t) {
  // y must be uniform 0-1
  let ret = 0;
  if (t === 1) {
    ret = 1;
  } else if (t === 2) {
    ret = Math.log(y / (1 - y));
    if (isNaN(ret)) {}
  } else if (t === 3) {
    ret = (y - 0.5) * Math.log(y / (1 - y));
    if (isNaN(ret)) {
      console.log("ret when t3 ", y, ret);
    }
  } else if (t === 4) {
    ret = y - 0.5;
    if (isNaN(ret)) {}
  } else if (t >= 5 && t % 2 === 1) {
    ret = Math.pow(y - 0.5, Math.floor((t - 1) / 2));
    if (isNaN(ret)) {}
  } else if (t >= 6 && t % 2 === 0) {
    ret = Math.pow(y - 0.5, Math.floor((t - 1) / 2)) * Math.log(y / (1 - y));
    if (isNaN(ret)) {}
  }
  return ret;
}

function flatten(arr) {
  return arr.reduce(function (flat, toFlatten) {
    return flat.concat(Array.isArray(toFlatten) ? flatten(toFlatten) : toFlatten);
  }, []);
}
