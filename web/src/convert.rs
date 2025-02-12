use bsxplorer_ci::convert::{ConvertReportType, IpcCompression};
use std::collections::HashMap;
use std::path::PathBuf;
use stylist::Style;
use web_sys::HtmlInputElement;
use yew::prelude::*;

#[derive(Clone, Debug, PartialEq, Properties)]
pub struct Config {
    #[prop_or_default]
    pub input_files: Vec<String>,
    #[prop_or_default]
    pub output_dir: String,
    #[prop_or(ConvertReportType::Bsx)]
    pub from_type: ConvertReportType,
    #[prop_or(ConvertReportType::Bsx)]
    pub into_type: ConvertReportType,
    #[prop_or(Some(IpcCompression::None))]
    pub ipc_compression: Option<IpcCompression>,
    #[prop_or(1)]
    pub threads: u32,
    #[prop_or(false)]
    pub low_memory: bool,
    #[prop_or(10_000)]
    pub chunk_size: u32,
    #[prop_or(None)]
    pub fasta_path: Option<String>,
    #[prop_or(None)]
    pub fai_path: Option<String>,
    #[prop_or(16)]
    pub batch_per_read: u32,
    #[prop_or(1_000_000)]
    pub batch_size: u32,
    #[prop_or(ViewState::Config)]
    pub state: ViewState,
    #[prop_or(0)]
    pub progress: u32,
}

#[derive(PartialEq, Eq, Debug, Clone)]
pub enum ViewState {
    Config,
    Confirm,
    Progress,
    Done,
}

/* Define Props with `UseStateHandle` */
#[derive(Properties, PartialEq)]
pub struct MyProps {
    pub value: UseStateHandle<u32>,
}

#[function_component(MyComponent)]
fn my_component(props: &MyProps) -> Html {
    let value_handle = props.value.clone();

    let on_increment = Callback::from(move |_| {
        let new_value = *value_handle + 1;
        value_handle.set(new_value);
    });

    html! {
        <div>
            <p>{ format!("Current value: {}", *props.value) }</p>
            <button onclick={on_increment}>{ "Increment" }</button>
        </div>
    }
}

#[derive(Properties, PartialEq)]
pub struct InputFilesProps {
    pub files: Vec<String>,
    pub invalid: HashMap<usize, String>,
    pub on_add_filepath: Callback<MouseEvent, ()>,
    pub on_remove_filepath: Callback<usize, ()>,
    pub on_change_filepath: Callback<(usize, InputEvent), ()>,
}

#[function_component(InputFilesForm)]
pub fn input_files_form(props: &InputFilesProps) -> Html {
    html! {
        <div class="field">
            // Field description
            <label for="input_files">
                { "Input Files" }
                <span class="tooltip-icon" data-tooltip="Enter file paths manually."> {"?"} </span>
            </label>
            // Files container
            <div class="file-paths">
                {
                    for props.files.iter().enumerate().map(|(index, path)| html! {
                        // File entry
                        <div class="file-path-entry">
                            // Input block
                            <input
                                type="text"
                                class={if props.invalid.get(&index).is_some() { "error" } else { "" }}
                                placeholder="Enter file path..."
                                value={path.clone()}
                                oninput={
                                    props.on_change_filepath.reform(move |e: InputEvent| (index, e))
                                }
                            />
                            // Remove button
                            <button
                                class="remove"
                                onclick={
                                    props.on_remove_filepath.reform(move |e: MouseEvent| {
                                        e.prevent_default();
                                        index
                                    })
                                }
                            > {"✖"}
                            </button>
                            // Error popup
                            if let Some(msg) = props.invalid.get(&index) {
                                <div class="error-popup">{ msg }</div>
                            }
                        </div>
                    })
                }
            </div>
            // Add file button
            <button
                class="add"
                onclick={props.on_add_filepath.clone()}
            > {"➕ Add File Path"}
            </button>
        </div>
    }
}

mod callbacks {
    use super::*;
    use gloo_timers::callback::Interval;
    use paste::paste;
    use web_sys::HtmlInputElement;

    pub fn on_change_file_path(config: UseStateHandle<Config>) -> Callback<(usize, InputEvent)> {
        let config = config.clone();
        Callback::from(move |(index, e): (usize, InputEvent)| {
            let input: HtmlInputElement = e.target_unchecked_into();
            let mut new_config = (*config).clone();
            new_config.input_files[index] = input.value();
            config.set(new_config);
        })
    }

    pub fn on_add_file_path(config: UseStateHandle<Config>) -> Callback<MouseEvent> {
        let config = config.clone();
        Callback::from(move |e: MouseEvent| {
            e.prevent_default();
            let mut new_config = (*config).clone();
            new_config.input_files.push("".to_string()); // Add empty entry
            config.set(new_config);
        })
    }

    pub fn on_remove_file_path(config: UseStateHandle<Config>) -> Callback<usize> {
        let config = config.clone();
        Callback::from(move |index: usize| {
            if config.input_files.len() > 1 {
                let mut new_config = (*config).clone();
                new_config.input_files.remove(index);
                config.set(new_config);
            }
        })
    }

    macro_rules! change_config {
        ($value:ident, $body:expr) => {
            paste! {
            pub fn [<on_change_ $value>](config: UseStateHandle<Config>) -> Callback<Event> {
                let config = config.clone();
                    Callback::from(move |e: Event| {
                        let res = $body(e);

                        let mut new_config = (*config).clone();
                        new_config.$value = res;
                        config.set(new_config);
                    })
                }
            }
        };
    }

    change_config!(output_dir, |e: Event| {
        let input: HtmlInputElement = e.target_unchecked_into();
        if let Ok(val) = input.value().parse::<String>() {
            val
        } else {
            "".to_string()
        }
    });

    change_config!(low_memory, |e: Event| {
        let input: HtmlInputElement = e.target_unchecked_into();
        input.checked()
    });

    change_config!(threads, |e: Event| {
        let input: HtmlInputElement = e.target_unchecked_into();
        if let Ok(num) = input.value().parse::<u32>() {
            num
        } else {
            1
        }
    });
    change_config!(chunk_size, |e: Event| {
        let input: HtmlInputElement = e.target_unchecked_into();
        if let Ok(num) = input.value().parse::<u32>() {
            num
        } else {
            10_000
        }
    });
    change_config!(batch_per_read, |e: Event| {
        let input: HtmlInputElement = e.target_unchecked_into();
        if let Ok(num) = input.value().parse::<u32>() {
            num
        } else {
            16
        }
    });
    change_config!(batch_size, |e: Event| {
        let input: HtmlInputElement = e.target_unchecked_into();
        if let Ok(num) = input.value().parse::<u32>() {
            num
        } else {
            1_000_000
        }
    });

    change_config!(into_type, |e: Event| {
        let select: web_sys::HtmlSelectElement = e.target_unchecked_into();
        let value = select.value();
        match value.as_str() {
            "Bsx" => ConvertReportType::Bsx,
            "Bismark" => ConvertReportType::Bismark,
            "CgMap" => ConvertReportType::CgMap,
            "BedGraph" => ConvertReportType::BedGraph,
            "Coverage" => ConvertReportType::Coverage,
            _ => ConvertReportType::Bsx,
        }
    });
    change_config!(from_type, |e: Event| {
        let select: web_sys::HtmlSelectElement = e.target_unchecked_into();
        let value = select.value();
        match value.as_str() {
            "Bsx" => ConvertReportType::Bsx,
            "Bismark" => ConvertReportType::Bismark,
            "CgMap" => ConvertReportType::CgMap,
            "BedGraph" => ConvertReportType::BedGraph,
            "Coverage" => ConvertReportType::Coverage,
            _ => ConvertReportType::Bsx,
        }
    });
    change_config!(ipc_compression, |e: Event| {
        let select: web_sys::HtmlSelectElement = e.target_unchecked_into();
        let value = select.value();
        Some(match value.as_str() {
            "LZ4" => IpcCompression::LZ4,
            "ZSTD" => IpcCompression::ZSTD,
            "None" => IpcCompression::None,
            _ => IpcCompression::None,
        })
    });

    change_config!(fasta_path, |e: Event| {
        let input: HtmlInputElement = e.target_unchecked_into();
        let files = input.files();
        let fasta = if let Some(files_list) = files {
            if files_list.length() > 0 {
                if let Some(file) = files_list.get(0) {
                    Some(file.name())
                } else {
                    None
                }
            } else {
                None
            }
        } else {
            None
        };
        fasta
    });

    change_config!(fai_path, |e: Event| {
        let input: HtmlInputElement = e.target_unchecked_into();
        let files = input.files();
        let fasta = if let Some(files_list) = files {
            if files_list.length() > 0 {
                if let Some(file) = files_list.get(0) {
                    Some(file.name())
                } else {
                    None
                }
            } else {
                None
            }
        } else {
            None
        };
        fasta
    });

    pub fn on_back(config: UseStateHandle<Config>) -> Callback<MouseEvent> {
        Callback::from(move |e: MouseEvent| {
            e.prevent_default();
            // In a complete app, you might validate & update `config` here.
            let mut new_config = (*config).clone();
            new_config.state = ViewState::Config;
            config.set(new_config);
        })
    }

    pub fn on_confirm(config: UseStateHandle<Config>) -> Callback<MouseEvent> {
        Callback::from(move |_| {
            let mut new_config = (*config).clone();
            new_config.state = ViewState::Progress;
            config.set(new_config);
            // Start an interval timer to update progress.
            let mut new_config = (*config).clone();

            let interval = Interval::new(100, move || {
                let current = new_config.progress;
                if current < 100 {
                    new_config.progress = current + 5;
                } else {
                    new_config.state = ViewState::Done;
                }
            });
            // Keep the interval alive.
            interval.forget();
        })
    }

    pub fn on_submit(config: UseStateHandle<Config>) -> Callback<SubmitEvent> {
        Callback::from(move |e: SubmitEvent| {
            e.prevent_default();
            // In a complete app, you might validate & update `config` here.
            let mut new_config = (*config).clone();
            new_config.state = ViewState::Confirm;
            config.set(new_config);
        })
    }
}

pub mod states {
    use super::*;
    use crate::convert::Config;
    use yew::html::IntoPropValue;
    use yew::{html, props, Html, UseStateHandle};

    #[function_component]
    pub fn ConfigView(config: &Config) -> Html {
        fn error_condition(condition: bool) -> &'static str {
            if condition {
                "error"
            } else {
                ""
            }
        }
        let is_confirm_disabled = config.input_files.is_empty()
            || config.output_dir.is_empty()
            || (matches!(
                config.from_type,
                ConvertReportType::BedGraph | ConvertReportType::Coverage
            ) && (config.fasta_path.is_none() || config.fai_path.is_none()))
            || (matches!(config.from_type, ConvertReportType::Bsx)
                && (config.fasta_path.is_none() && config.fai_path.is_none()));
        let fasta_error = (matches!(
            config.from_type,
            ConvertReportType::BedGraph | ConvertReportType::Coverage
        ) && (config.fasta_path.is_none() || config.fai_path.is_none()))
            || (matches!(config.from_type, ConvertReportType::Bsx)
                && (config.fasta_path.is_none() && config.fai_path.is_none()));

        let files = Vec::<String>::new();
        let invalid = HashMap::<usize, String>::new();
        let files_handle = use_state(|| files);
        let invalid_handle = use_state(|| invalid);

        let on_add_filepath = {
            let files_handle = files_handle.clone();
            Callback::from(move |e: MouseEvent| {
                e.prevent_default();
                let mut new_files = (*files_handle).clone();
                new_files.push("".to_string());
                files_handle.set(new_files);
            })
        };

        let on_remove_filepath = {
            let files_handle = files_handle.clone();
            Callback::from(move |index: usize| {
                let mut new_files = (*files_handle).clone();
                new_files.remove(index);
                files_handle.set(new_files);
            })
        };

        let on_change_filepath = {
            let files_handle = files_handle.clone();
            let invalid_handle = invalid_handle.clone();

            Callback::from(move |(index, e): (usize, InputEvent)| {
                let input: HtmlInputElement = e.target_unchecked_into();
                let value = input.value();

                if value.is_empty() {
                    // Empty string
                    let mut new_invalid = (*invalid_handle).clone();
                    let entry = new_invalid.entry(index).or_default();
                    entry.clear();
                    entry.push_str("Input file path can not be empty");
                } else {
                    let path = PathBuf::from(value);
                    if !path.exists() {
                        web_sys::console::log_1(&format!("Updated file list: {:?}", path).into());
                        let mut new_invalid = (*invalid_handle).clone();
                        let entry = new_invalid.entry(index).or_default();
                        entry.clear();
                        entry.push_str("Input file does not exist");
                    // Not file
                    } else if !path.is_file() {
                        let mut new_invalid = (*invalid_handle).clone();
                        let entry = new_invalid.entry(index).or_default();
                        entry.clear();
                        entry.push_str("Input is not a file");
                    } else {
                        // ACTUALLY set the file path
                        let mut new_files = (*files_handle).clone();
                        new_files[index] = input.value();
                        files_handle.set(new_files);
                    }
                }
            })
        };

        html! {
            <form class="container"
            onsubmit={|_| println!("1") }
            >
            // --- Input files (multiple) ---
            <div class="field">
            <label for="input_files">
                { "Input Files" }
                <span class="tooltip-icon" data-tooltip="Enter file paths manually."> {"?"} </span>
            </label>

                <InputFilesForm {on_add_filepath} {on_change_filepath} {on_remove_filepath} files={(*files_handle).clone()} invalid={(*invalid_handle).clone()}/>

            </div>

            <button disabled={is_confirm_disabled} type="submit" class="submit">{"Submit"}</button>
        </form>
        }
    }

    // html! {
    //         <form class="container" onsubmit={callbacks::on_submit(config.clone())}>
    //         // --- Input files (multiple) ---
    //         <div class="field">
    //         <label for="input_files">
    //             { "Input Files" }
    //             <span class="tooltip-icon" data-tooltip="Enter file paths manually."> {"?"} </span>
    //         </label>
    //
    //             <div class="file-paths">
    //                 {
    //                     for config.input_files.iter().enumerate().map(|(index, path)| html! {
    //                         <div class="file-path-entry">
    //                             <input
    //                                 type="text"
    //                                 placeholder="Enter file path..."
    //                                 value={path.clone()}
    //                                 oninput={callbacks::on_change_file_path(config.clone()).reform(move |e: InputEvent| (index, e))}
    //                             />
    //                             <button class="remove" onclick={callbacks::on_remove_file_path(config.clone()).reform(move |e: MouseEvent| {
    //                                 e.prevent_default();
    //                                 index
    //                             })}>{"✖"}</button>
    //                         </div>
    //                     })
    //                 }
    //             </div>
    //
    //             <button class="add" onclick={callbacks::on_add_file_path(config.clone())}>{"➕ Add File Path"}</button>
    //         </div>
    //
    //         // --- Output directory selection ---
    //         <div class="field">
    //             <label for="output_dir">
    //                 {"Output prefix"}
    //                 <span class="tooltip-icon" data-tooltip="Prefix for the generated output files."> {"?"} </span>
    //             </label>
    //             <input
    //                 id="output_dir"
    //                 type="text"
    //                 onchange={callbacks::on_change_output_dir(config.clone())}
    //                 class={error_condition(config.output_dir.len() == 0)}
    //             />
    //         </div>
    //
    //         // --- From Type selection ---
    //         <div class="field">
    //             <label for="from_type">{"From Type"}</label>
    //             <select id="from_type" onchange={callbacks::on_change_from_type(config.clone())}>
    //                 <option value="Bsx" selected={config.from_type == ConvertReportType::Bsx}>{"Bsx"}</option>
    //                 <option value="Bismark" selected={config.from_type == ConvertReportType::Bismark}>{"Bismark"}</option>
    //                 <option value="CgMap" selected={config.from_type == ConvertReportType::CgMap}>{"CgMap"}</option>
    //                 <option value="BedGraph" selected={config.from_type == ConvertReportType::BedGraph}>{"BedGraph"}</option>
    //                 <option value="Coverage" selected={config.from_type == ConvertReportType::Coverage}>{"Coverage"}</option>
    //             </select>
    //         </div>
    //
    //         // --- Into Type selection ---
    //         <div class="field">
    //             <label for="into_type">{"Into Type"}</label>
    //             <select id="into_type" onchange={callbacks::on_change_into_type(config.clone())} >
    //                 <option value="Bsx" selected={config.into_type == ConvertReportType::Bsx}>{"Bsx"}</option>
    //                 <option value="Bismark" selected={config.into_type == ConvertReportType::Bismark}>{"Bismark"}</option>
    //                 <option value="CgMap" selected={config.into_type == ConvertReportType::CgMap}>{"CgMap"}</option>
    //                 <option value="BedGraph" selected={config.into_type == ConvertReportType::BedGraph}>{"BedGraph"}</option>
    //                 <option value="Coverage" selected={config.into_type == ConvertReportType::Coverage}>{"Coverage"}</option>
    //             </select>
    //         </div>
    //
    //         // --- IPC Compression (only shown if into_type == Bsx) ---
    //         { if let ConvertReportType::Bsx = config.into_type {
    //             html! {
    //                 <div class="field">
    //                     <label for="ipc_compression">{"IPC Compression"}</label>
    //                     <select id="ipc_compression" onchange={callbacks::on_change_ipc_compression(config.clone())}>
    //                         <option value="LZ4">{"LZ4"}</option>
    //                         <option value="ZSTD">{"ZSTD"}</option>
    //                         <option value="None" selected={matches!(config.ipc_compression, Some(IpcCompression::None))}>{"None"}</option>
    //                     </select>
    //                 </div>
    //             }
    //         } else {
    //             html! {}
    //         }}
    //
    //         // --- Threads input ---
    //         <div class="field">
    //             <label for="threads">
    //                 {"Threads"}
    //                 <span class="tooltip-icon" data-tooltip="Number of threads to use."> {"?"} </span>
    //             </label>
    //             <input
    //                 id="threads"
    //                 type="number"
    //                 min="1"
    //                 max="16"
    //                 value={config.threads.to_string()}
    //                 onchange={callbacks::on_change_threads(config.clone())}
    //             />
    //         </div>
    //
    //         // --- Low Memory checkbox ---
    //         <div class="field checkbox-field">
    //             <input
    //                 id="low_memory"
    //                 type="checkbox"
    //                 checked={config.low_memory}
    //                 onchange={callbacks::on_change_low_memory(config.clone())}
    //             />
    //             <label for="low_memory">
    //                 {"Low Memory"}
    //                 <span class="tooltip-icon" data-tooltip="Use less RAM, but elongate computation."> {"?"} </span>
    //             </label>
    //         </div>
    //
    //         // --- Chunk Size input ---
    //         <div class="field">
    //             <label for="chunk_size">
    //                 {"Chunk Size"}
    //                 <span class="tooltip-icon" data-tooltip="Number of rows in the output batches (Important when converting to bsx format)."> {"?"} </span>
    //             </label>
    //             <input
    //                 id="chunk_size"
    //                 type="number"
    //                 value={config.chunk_size.to_string()}
    //                 onchange={callbacks::on_change_chunk_size(config.clone())}
    //             />
    //         </div>
    //
    //         // --- FASTA Path: required if from_type is BedGraph or Coverage ---
    //         { if matches!(config.from_type, ConvertReportType::BedGraph | ConvertReportType::Coverage) {
    //             html! {
    //                 <div class="field">
    //                     <label for="fasta_path">
    //                         <span class="tooltip-icon" data-tooltip="Path to the reference sequence file. Obligatory when converting BedGraph or Coverage."> {"?"} </span>
    //                         {"FASTA Path (required)"}
    //                     </label>
    //                     <input id="fasta_path" type="file" class={error_condition(fasta_error)} onchange={callbacks::on_change_fasta_path(config.clone())} />
    //                 </div>
    //             }
    //         } else {
    //             html! {
    //                 <div class="field">
    //                     <label for="fasta_path">
    //                         {"FASTA Path (optional)"}
    //                         <span class="tooltip-icon" data-tooltip="Path to the reference sequence file. Obligatory when converting BedGraph or Coverage."> {"?"} </span>
    //                     </label>
    //                     <input id="fasta_path" type="file" onchange={callbacks::on_change_fasta_path(config.clone())} />
    //                 </div>
    //             }
    //         }}
    //
    //         // --- FAI Path: required if from_type is BedGraph/Coverage or into_type is Bsx ---
    //         { if matches!(config.from_type, ConvertReportType::BedGraph | ConvertReportType::Coverage) || matches!(config.into_type, ConvertReportType::Bsx) {
    //             html! {
    //                 <div class="field" >
    //                     <label for="fai_path">
    //                         {"FAI Path (required)"}
    //                         <span class="tooltip-icon" data-tooltip="Path to the fasta index. Obligatory when converting BedGraph or Coverage."> {"?"} </span>
    //                     </label>
    //                     <input id="fai_path" class={error_condition(fasta_error)}  type="file" onchange={callbacks::on_change_fai_path(config.clone())} />
    //                 </div>
    //             }
    //         } else {
    //             html! {
    //                 <div class="field">
    //                     <label for="fai_path">
    //                         {"FAI Path (optional)"}
    //                         <span class="tooltip-icon" data-tooltip="Path to the fasta index. Obligatory when converting BedGraph or Coverage."> {"?"} </span>
    //                     </label>
    //                     <input id="fai_path" type="file" onchange={callbacks::on_change_fai_path(config.clone())} />
    //                 </div>
    //             }
    //         }}
    //
    //         // --- Batch Per Read input ---
    //         <div class="field">
    //             <label for="batch_per_read">
    //                 {"Batch Per Read"}
    //                 <span class="tooltip-icon" data-tooltip="Number of batches to read simultaneously. Affects RAM usage."> {"?"} </span>
    //             </label>
    //             <input
    //                 id="batch_per_read"
    //                 type="number"
    //                 value={config.batch_per_read.to_string()}
    //                 onchange={callbacks::on_change_batch_per_read(config.clone())}
    //             />
    //         </div>
    //
    //         // --- Batch Size input ---
    //         <div class="field">
    //             <label for="batch_size">
    //                 {"Batch Size"}
    //                 <span class="tooltip-icon" data-tooltip="Read batch size. Affects RAM usage."> {"?"} </span>
    //             </label>
    //             <input
    //                 id="batch_size"
    //                 type="number"
    //                 value={config.batch_size.to_string()}
    //                 onchange={callbacks::on_change_batch_size(config.clone())}
    //             />
    //         </div>
    //
    //         <button disabled={is_confirm_disabled} type="submit" class="submit">{"Submit"}</button>
    //     </form>
    //     }

    pub fn confirm(config: UseStateHandle<Config>) -> Html {
        html! {
            <div class="container">
                <h2>{ "Confirm Your Configuration" }</h2>
                <label for="item-list">{"Input Files"}</label>
                <ul class="item-list">
                        { for config.input_files.clone().into_iter().map(|v| html! {<li> { v } </li>}) }
                </ul>
                <p>{ format!("Output prefix: {}", config.output_dir) }</p>
                <p>{ format!("From Type: {:?}", config.from_type) }</p>
                <p>{ format!("Into Type: {:?}", config.into_type) }</p>
                if { matches!(config.into_type, ConvertReportType::Bsx) } {
                    <p> {format!("IPC file compression: {:?}", config.into_type)} </p>
                }
                <p>{ format!("Threads: {}", config.threads) }</p>
                <p>{ format!("Low-memory: {}", config.low_memory) }</p>
                <p>{ format!("Chunk size: {}", config.chunk_size) }</p>
                if { config.fasta_path.is_some() } {
                   <p> {format!("Fasta path: {}", config.fasta_path.as_ref().unwrap())} </p>
                }
                if { config.fai_path.is_some() } {
                   <p> {format!("Fasta path: {}", config.fai_path.as_ref().unwrap())} </p>
                }
                <p>{ format!("Batch per read: {}", config.batch_per_read) }</p>
                <p>{ format!("Read batch size: {}", config.batch_size) }</p>
                <div class="button-container">
                    <button class="deny" onclick={callbacks::on_back(config.clone())}>{ "Deny" }</button>
                    <button class="confirm" onclick={callbacks::on_back(config.clone())}>{ "Confirm" }</button>
                </div>
            </div>
        }
    }

    pub fn progress(config: UseStateHandle<Config>) -> Html {
        html! {
            <div class="container">
                <h2>{ "Processing..." }</h2>
                <progress max="100" value={config.progress.to_string()}>{ format!("{}%", config.progress) }</progress>
            </div>
        }
    }

    pub fn done(config: UseStateHandle<Config>) -> Html {
        html! {
            <div class="container">
                <h2>{ "Process Complete" }</h2>
                <button onclick={callbacks::on_back(config).clone()}>{ "Back to Config" }</button>
            </div>
        }
    }
}

pub fn style() -> Style {
    stylist::style!(
        r#"
        
        "#
    )
    .unwrap()
}
