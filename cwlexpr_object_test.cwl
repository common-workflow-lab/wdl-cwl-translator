cwlVersion: v1.0
class: ExpressionTool
requirements:
 InlineJavascriptRequirement:
   expressionLib: [ $include: map.js ]

inputs: []

expression: |
  ${
    var samOutputNames = new Map([ ["BAM SortedByCoordinate", "sortedByCoord.out.bam"] ]);
    return { "only": "prefix" + "Aligned." + samOutputNames.get("BAM SortedByCoordinate")} ;
   }

# Above uses Map from ECMAScript 6 (ECMAScript 2015)

# Below uses an Object object, compatible with the 5th Edition of ECMAScript

#  ${
#    var samOutputNames = {"BAM SortedByCoordinate": "sortedByCoord.out.bam"};
#    return {"only": "prefix" + "Aligned." + samOutputNames["BAM SortedByCoordinate"]};
#   }

outputs:
  only: string

