
# Converts severus vcf to have SVTYPE=DEL instead of BND

import pysam

vcf = pysam.VariantFile('-')
vcf_out = pysam.VariantFile('-', 'w', header=vcf.header)
for rec in vcf:
    if 'SVTYPE' not in rec.info or rec.info['SVTYPE'] != 'BND':
        vcf_out.write(rec)
        continue
    if 'DETAILED_TYPE' in rec.info and rec.info['DETAILED_TYPE'] != 'DEL':
        vcf_out.write(rec)
        continue
    if 'MATE_ID' not in rec.info:
        vcf_out.write(rec)
        continue
    if rec.info['MATE_ID'].endswith('_2'):
        continue
    
    rec.info['SVTYPE'] = 'DEL'

    pos1 = rec.pos
    pos2 = rec.alts[0].replace('N', '').replace(']', '').replace('[', '').split(':')[1]
    pos1, pos2 = sorted([pos1, int(pos2)])
    
    rec.pos = pos1
    rec.alts = ('<DEL>',)

    vcf_out.write(rec)


vcf_out.close()
vcf.close()

