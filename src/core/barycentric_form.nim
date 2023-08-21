import 
 constantine/math/config/[type_ff, curves],
 constantine/math/elliptic/ec_twistededwards_projective,
 constantine/math/arithmetic/[finite_fields, bigints,bigints_montgomery]


type 
 PrecomputedWeights = object
  barycentricWeights: seq[ECP_TwEdwards_Prj[Fr[Bandersnatch]]]
  invertedDomain: seq[ECP_TwEdwards_Prj[Fr[Bandersnatch]]]

const
 DOMAIN: uint64 = 256

# proc UpdatePrecomputedWeights(newWeight : ECP_TwEdwards_Prj[Fr[Bandersnatch]], newDomain : ECP_TwEdwards_Prj[Fr[Bandersnatch]]):
 
var PrecomputedWeightsObj: PrecomputedWeights

proc ComputeBaryCentricWeights (element : uint64) : ECP_TwEdwards_Prj[Fr[Bandersnatch]] = 
 if element > DOMAIN:
  echo"The domain is [0,255], and $element is not in the domain"

 var domain_element_fr: ECP_TwEdwards_Prj[Fr[Bandersnatch]]
#  var conv_domain_element_fr: uint64 = cast[uint64](domain_element_fr)

 var total : ECP_TwEdwards_Prj[Fr[Bandersnatch]]

 total.x.setOne()
 total.y.setOne()
 total.z.setOne()

 for i in uint64(0)..DOMAIN:
  if i == element:
    continue

  var i_fr: ECP_TwEdwards_Prj[Fr[Bandersnatch]] = cast[ECP_TwEdwards_Prj[Fr[Bandersnatch]]](i)
  # var conv_i_fr: uint64 = cast[uint64](i_fr)
  var temp : ECP_TwEdwards_Prj[Fr[Bandersnatch]]
  temp.diff(domain_element_fr,i_fr)

  total.x.prod(total.x,temp.x)
  total.y.prod(total.y,temp.y)
  total.z.prod(total.z,temp.z)
  
 return total



proc NewPreComputedWeights() : PrecomputedWeights =
 var midpoint: uint64 = DOMAIN
 var barycentricWeightsInst: seq[ECP_TwEdwards_Prj[Fr[Bandersnatch]]] = newSeq[ECP_TwEdwards_Prj[Fr[Bandersnatch]]](midpoint * 2)
 
 for i in uint64(0)..midpoint:
  var weights : ECP_TwEdwards_Prj[Fr[Bandersnatch]] = ComputeBaryCentricWeights(i)

  var inverseWeights : ECP_TwEdwards_Prj[Fr[Bandersnatch]]

  inverseWeights.x.inv(weights.x)
  inverseWeights.y.inv(weights.y)
  inverseWeights.z.inv(weights.z)


  barycentricWeightsInst[i] = weights
  barycentricWeightsInst[i+midpoint] = inverseWeights


  midpoint = DOMAIN - 1
  var invertedDomain: seq[ECP_TwEdwards_Prj[Fr[Bandersnatch]]] = newSeq[ECP_TwEdwards_Prj[Fr[Bandersnatch]]](midpoint * 2)

  for i in uint64(0)..DOMAIN:
   var k: ECP_TwEdwards_Prj[Fr[Bandersnatch]] = cast[ECP_TwEdwards_Prj[Fr[Bandersnatch]]](i)

   k.x.inv(k.x)
   k.y.inv(k.y)
   k.z.inv(k.z)

   var neg_k : ECP_TwEdwards_Prj[Fr[Bandersnatch]]

   var zero : ECP_TwEdwards_Prj[Fr[Bandersnatch]]

   zero.x.setZero()
   zero.y.setZero()
   zero.z.setZero()

   neg_k.diff(zero, k)

   invertedDomain[i-1] = k
   invertedDomain[(i-1) + midpoint] = neg_k
   
  
   PrecomputedWeightsObj.barycentricWeights = barycentricWeightsInst
   PrecomputedWeightsObj.invertedDomain = invertedDomain

   return PrecomputedWeightsObj

proc BatchInversion(points : seq[ECP_TwEdwards_Prj[Fr[Bandersnatch]]]) : seq[ECP_TwEdwards_Prj[Fr[Bandersnatch]]] =
 var result : array[len(points),ECP_TwEdwards_Prj[Fr[Bandersnatch]]]

proc ComputeBaryCentricCoefficients( point : ECP_TwEdwards_Prj[Fr[Bandersnatch]]): seq[ECP_TwEdwards_Prj[Fr[Bandersnatch]]] =
 var lagrangeEval : array[DOMAIN, ECP_TwEdwards_Prj[Fr[Bandersnatch]]]

 for i in uint64(0)..DOMAIN:
  var weight = PrecomputedWeightsObj.barycentricWeights[i]
  var i_fr: ECP_TwEdwards_Prj[Fr[Bandersnatch]] = cast[ECP_TwEdwards_Prj[Fr[Bandersnatch]]](i)
  
  lagrangeEval[i].diff(point, i_fr)

  lagrangeEval[i].x.prod(lagrangeEval[i].x,weight.x)
  lagrangeEval[i].y.prod(lagrangeEval[i].y,weight.y)
  lagrangeEval[i].z.prod(lagrangeEval[i].z,weight.z)
  
var totalProd : ECP_TwEdwards_Prj[Fr[Bandersnatch]]

totalProd.x.setOne()
totalProd.y.setOne()
totalProd.z.setOne()




   