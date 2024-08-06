# A. Irles
# IFIC - CSIC/UV
# https://aitanatop.ific.uv.es/aitanatop/luxe/
# https://aitanatop.ific.uv.es/aitanatop/detector-rd/
# August 2024

# SimCalo2CaloHit

Processor that will read SimCalorimeterHits and convert them to CalorimeterHits in MIP units, using the input from MIPCalibration.

This processor should also produce a digitized signal - in first approach we will simply smear the energy by the S/N value of the ADC for MIPs, expected to be 10.


Work-in-progress