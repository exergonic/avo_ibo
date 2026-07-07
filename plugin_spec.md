# Plugin Info

Taken from [How Plugins Are Run](https://avogadro.cc/develop/plugins/args.html)

## How plugins are run

- Plugin functionality is invoked by Avogadro via two different APIs:
  1. Script API
  2. Package API
- Each feature that a plugin offers has an `identifier` as listed in the plugin’s metadata file
- Features of a script plugin are invoked by Avogadro using `pixi run python plugins/python/<plugin>/<identifier>.py [OPTIONS]`
- Features of packages (including Pixi projects) are invoked using `pixi run avogadro-<plugin> <identifier> [OPTIONS]`
  - e.g. for the `mass` feature of a plugin called `demo` the full command would be `pixi run avogadro-demo mass`
  - The plugin defines an entry point in their `pyproject.toml` for the plugin’s name prefixed with `avogadro-` e.g. for a `demo` plugin:

```toml
[project.scripts]
avogadro-demo = "demo:main"
```

- A plugin is required to correctly process certain option flags (by which it is meant that the plugin should not crash or error):
  - `--lang` – though actually changing behaviour based on the accompanying localization string is not obligatory
  - `--debug` – again, a change in behaviour is not obligatory, though plugins are encouraged to provide more detailed logs and error messages if this flag is passed
  - `--print-options` – *only* if the plugin indicates that the interface of a feature’s dialog is to be specified *dynamically* using `user-options = "dynamic"`, it must respond appropriately to the flag, otherwise it does not need to as it will never be passed
  - `electrostatic-models` features must respond appropriately (according to what they support) to:
    - `--charges`
    - `--potential`
  - `file-formats` features must respond appropriately (according to what they support) to:
    - `--read`
    - `--write`
- All other option flags used in the plugin API until now are deprecated
  - `--run-command` and similar are no longer used; if a feature type only ever responds to one kind of request, no flag is needed, it is assumed to be the desired request
  - The information provided by things like `--display-name` and `--menu-path`, is now contained in the static metadata
- Avogadro will now *always* pass a JSON object as the main argument (following the feature identifier) to all feature types – for details of the contents, see below
  - This is no longer passed on `stdin`; streams of subsequent data are, however, if relevant for the plugin type
- Avogadro will now *always* accept an JSON object as the output of the plugin to `stdout` – for details of the contents, see below
  - If Avogadro expects a stream of output on `stdout`, it will first expect a JSON object as the initial output before the stream, and will handle any subsequent JSON object on `stdout` as an indication that the plugin process should be stopped

> NOTE: several examples locally in `C:\Users\mccan\AppData\Local\OpenChemistry\Avogadro\plugins\` 
>
> the `avo_xtb` plugin also returns orbital much in the same way we wish to return orbitals in our plugin, when the `orbitals` option is invoked.
>
> we may wish to compare how `avo_xtb` creates and returns its `molden_string`

## The input JSON

| Key                                                   | JSON data type      | Value                                                       | Feature types                       | Applicable                                       |
| ----------------------------------------------------- | ------------------- | ----------------------------------------------------------- | ----------------------------------- | ------------------------------------------------ |
| `config`                                              | `object`            | The user’s stored configuration of the plugin               | All                                 | `user-config` is specified                       |
| `cjson`                                               | `object`            | The current CJSON file in full                              | All except `file-formats`           | `input-cjson = true` or `input-format = "cjson"` |
| `xyz`, `cml`, `pdb`, `cif`, `sdf`, `mdl`, `mol`, etc. | `string`            | The contents of a file in the respective format             | All except `file-formats`           | `input-format = <fmt>`                           |
| `file`                                                | `string`            | A path to a file in the filesystem                          | `file-formats`                      | Always                                           |
| `options`                                             | `object`            | Key/value pairs for each option requested in `user-options` | `menu-commands`, `input-generators` | Always                                           |
| `points`                                              | `array` of `number` | The coordinates of the requested points                     | `electrostatic-models`              | Called with `--potentials`                       |
| `selectedAtoms`                                       | `array` of `number` | The indices of the atoms that are currently selected        | `menu-commands`                     | Always                                           |

## The output JSON

| Key                                                   | JSON data type      | Value                                                        | Feature types          | Applicable                 |
| ----------------------------------------------------- | ------------------- | ------------------------------------------------------------ | ---------------------- | -------------------------- |
| `message`                                             | `string`            | A message to display to the user in a pop-up dialog          | All                    | Always                     |
| `error`                                               | `string`            | Similar, but displays an error, and the plugin is stopped    | All                    | Always                     |
| `config`                                              | `object`            | The user’s configuration of the plugin as it should be saved | All                    | `user-config` is specified |
| `cjson`                                               | `object`            | A CJSON file to be loaded                                    | All?                   | `output-format = "cjson"`  |
| `xyz`, `cml`, `pdb`, `cif`, `sdf`, `mdl`, `mol`, etc. | `string`            | The contents of a file in the respective format              | All?                   | `output-format = <fmt>`    |
| `charges`                                             | `array` of `number` | The partial charges on each atom                             | `electrostatic-models` | Called with `--charges`    |
| `potentials`                                          | `array` of `number` | The electrostatic potential at each requested point          | `electrostatic-models` | Called with `--potentials` |
| `bond`                                                | `bool`              | Whether Avogadro should run automatic bond perception        | `menu-commands`        | Always                     |
| `selectedAtoms`                                       | `array` of `number` | The indices of the atoms that Avogadro should select         | `menu-commands`        | Always                     |
| `warnings`                                            | `array` of `string` | Non-fatal warnings to be shown in the generator’s dialog     | `input-generators`     | Always                     |
| `generatedFiles`                                      | `object`            | Details of several files                                     | `input-generators`     | Always                     |
| `mainFile`                                            | `string`            | The primary input file for a calculation                     | `input-generators`     | Always                     |
