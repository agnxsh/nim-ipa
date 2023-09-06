import 
 ../constantine/constantine/math/config/[type_ff, curves],
 ../constantine/constantine/math/elliptic/ec_twistededwards_projective,
 ../constantine/constantine/math/arithmetic/[finite_fields, bigints,bigints_montgomery]


type 
 PrecomputedWeights = object
  barycentricWeights: seq[ECP_TwEdwards_Prj[Fp[Bandersnatch]]]
  invertedDomain: seq[ECP_TwEdwards_Prj[Fp[Bandersnatch]]]

const
 DOMAIN: uint64 = 256

# proc UpdatePrecomputedWeights(newWeight : ECP_TwEdwards_Prj[Fp[Bandersnatch]], newDomain : ECP_TwEdwards_Prj[Fp[Bandersnatch]]):
 
var PrecomputedWeightsObj: PrecomputedWeights

proc ComputeBaryCentricWeights (element : uint64) : ECP_TwEdwards_Prj[Fp[Bandersnatch]] = 
 if element > DOMAIN:
  echo"The domain is [0,255], and $element is not in the domain"

 var domain_element_Fp: ECP_TwEdwards_Prj[Fp[Bandersnatch]]
#  var conv_domain_element_Fp: uint64 = cast[uint64](domain_element_Fp)

 var total : ECP_TwEdwards_Prj[Fp[Bandersnatch]]

 total.x.setOne()
 total.y.setOne()
 total.z.setOne()

 for i in uint64(0)..DOMAIN:
  if i == element:
    continue

  var i_Fp: ECP_TwEdwards_Prj[Fp[Bandersnatch]] = cast[ECP_TwEdwards_Prj[Fp[Bandersnatch]]](i)
  # var conv_i_Fp: uint64 = cast[uint64](i_Fp)
  var temp : ECP_TwEdwards_Prj[Fp[Bandersnatch]]
  temp.diff(domain_element_Fp,i_Fp)

  total.x.prod(total.x,temp.x)
  total.y.prod(total.y,temp.y)
  total.z.prod(total.z,temp.z)
  
 return total



proc NewPreComputedWeights() : PrecomputedWeights =
 var midpoint: uint64 = DOMAIN
 var barycentricWeightsInst {.noInit.} : seq[ECP_TwEdwards_Prj[Fp[Bandersnatch]]] = newSeq[ECP_TwEdwards_Prj[Fp[Bandersnatch]]](midpoint * 2)
 
 for i in uint64(0)..midpoint:
  var weights : ECP_TwEdwards_Prj[Fp[Bandersnatch]] = ComputeBaryCentricWeights(i)

  var inverseWeights : ECP_TwEdwards_Prj[Fp[Bandersnatch]]

  inverseWeights.x.inv(weights.x)
  inverseWeights.y.inv(weights.y)
  inverseWeights.z.inv(weights.z)


  barycentricWeightsInst[i].x = weights.x
  barycentricWeightsInst[i].y = weights.y
  barycentricWeightsInst[i].z = weights.z

  barycentricWeightsInst[i+midpoint].x = inverseWeights.x
  barycentricWeightsInst[i+midpoint].y = inverseWeights.y
  barycentricWeightsInst[i+midpoint].z = inverseWeights.z
  



  midpoint = DOMAIN - 1
  var invertedDomain: seq[ECP_TwEdwards_Prj[Fp[Bandersnatch]]] = newSeq[ECP_TwEdwards_Prj[Fp[Bandersnatch]]](midpoint * 2)

  for i in uint64(0)..DOMAIN:
   var k: ECP_TwEdwards_Prj[Fp[Bandersnatch]] = cast[ECP_TwEdwards_Prj[Fp[Bandersnatch]]](i)

   k.x.inv(k.x)
   k.y.inv(k.y)
   k.z.inv(k.z)

   var neg_k : ECP_TwEdwards_Prj[Fp[Bandersnatch]]

   var zero : ECP_TwEdwards_Prj[Fp[Bandersnatch]]

   zero.x.setZero()
   zero.y.setZero()
   zero.z.setZero()

   neg_k.diff(zero, k)

   invertedDomain[i-1].x = k.x
   invertedDomain[i-1].y = k.y
   invertedDomain[i-1].z = k.z

   invertedDomain[(i-1) + midpoint].x = neg_k.x
   invertedDomain[(i-1) + midpoint].y = neg_k.y
   invertedDomain[(i-1) + midpoint].z = neg_k.z
   
  
   PrecomputedWeightsObj.barycentricWeights = barycentricWeightsInst
   PrecomputedWeightsObj.invertedDomain = invertedDomain

   return PrecomputedWeightsObj

# proc BatchInversion(points : seq[ECP_TwEdwards_Prj[Fp[Bandersnatch]]]) : seq[ECP_TwEdwards_Prj[Fp[Bandersnatch]]] =
#  var result : array[len(points),ECP_TwEdwards_Prj[Fp[Bandersnatch]]]

func ComputeBaryCentricCoefficients* [PrecomputedWeightsObj]( point : ECP_TwEdwards_Prj[Fp[Bandersnatch]]): seq[ECP_TwEdwards_Prj[Fp[Bandersnatch]]] =
 var lagrangeEval : array[DOMAIN, ECP_TwEdwards_Prj[Fp[Bandersnatch]]]

 for i in uint64(0)..DOMAIN:
  var weight = PrecomputedWeightsObj.barycentricWeights[i]
  var i_Fp: ECP_TwEdwards_Prj[Fp[Bandersnatch]] = cast[ECP_TwEdwards_Prj[Fp[Bandersnatch]]](i)
  
  lagrangeEval[i].diff(point, i_Fp)

  lagrangeEval[i].x.prod(lagrangeEval[i].x,weight.x)
  lagrangeEval[i].y.prod(lagrangeEval[i].y,weight.y)
  lagrangeEval[i].z.prod(lagrangeEval[i].z,weight.z)
  
var totalProd : ECP_TwEdwards_Prj[Fp[Bandersnatch]]

totalProd.x.setOne()
totalProd.y.setOne()
totalProd.z.setOne()

for i in uint64(0)..DOMAIN:
 var i_fr {.noInit.} : ECP_TwEdwards_Prj[Fp[Bandersnatch]]




   