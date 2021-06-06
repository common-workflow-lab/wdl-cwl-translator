#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
import cwl_utils.parser_v1_2 as cwl

from ruamel import yaml

def main() -> None:
    """Generate a CWL object to match "cat-tool.cwl"."""
    inputs = []
    inputs.append(cwl.CommandInputParameter(id="input_bam", type="File"))
    inputs.append(cwl.CommandInputParameter(id="metrics_filename", type="string"))
    '''inputs = [
        cwl.CommandInputParameter(id="input_bam", type="File"),
        cwl.CommandInputParameter(id="metrics_filename", type="string")
        ]'''
    
    docker_requirement = [
        cwl.DockerRequirement(
            dockerPull="us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8",
        )
    ]

    hints = [
        cwl.ResourceRequirement(
            ramMin=3584,
        )
    ]

    outputs = [
        cwl.CommandOutputParameter(
            id="quality_yield_metrics",
            type="File",
            outputBinding=cwl.CommandOutputBinding(glob="$(inputs.metrics_filename)"),
        )
    ]

    arguments = [
        cwl.CommandLineBinding(
            valueFrom="INPUT=$(inputs.input_bam.path)"
        ),
        cwl.CommandLineBinding(
            valueFrom="OQ=true"
        ),
        cwl.CommandLineBinding(
            valueFrom="OUTPUT=$(inputs.metrics_filename)"
        )
    ]
    cat_tool = cwl.CommandLineTool(
        id="CollectQualityYieldMetrics",
        inputs=inputs,
        requirements=docker_requirement,
        hints=hints,
        outputs=outputs,
        cwlVersion="v1.0",
        baseCommand=['java', '-Xms2000m', '-jar', '/usr/picard/picard.jar', 'CollectQualityYieldMetrics' ],
        arguments=arguments,
    )

    print(yaml.main.round_trip_dump(cat_tool.save()))


if __name__ == "__main__":
    main()


'''python create_cwl_from_objects.py > result.cwl
cwltool --validate result.cwl'''